
pub mod common {

    use std::fmt;
    use chrono::{DateTime, Utc};
    use super::ascii::*;


    #[derive(Debug)]
    pub struct Reception {
        pub date_time: DateTime<Utc>,
        pub station: String,
        pub channel: usize,
        pub frequency: f64,
        pub level: isize,
        // avlc_block: Option
        pub acars_block: Option<Block>,
    }

    #[allow(dead_code)]
    #[derive(Debug)]
    pub struct Block {
        // B.PR-MAQ.H10.D11AJJ3242#DFBA47/A31947,1,1/AMDAR,REPORT/C1100,PR-MAQ,2337S,04639W/C203,102229,02370,0
        // ---------------------------------------
        // 0123456789012345678901234567890123456789
        // ||      || |||   |Flight Identifier (Only mandatory in Downlinks)
        // ||      || |||Message Sequence Number (Only mandatory in Downlinks)
        // ||      || |||Text (Optional)
        // ||      || ||Start of Text (Optional)
        // ||      || |Block Identifier
        // ||      ||Label
        // ||      |Technical Ack
        // ||Address
        // |Mode character
        // ----------------------------------------
        // AEEC SPECIFICATION 620-9 (PUBLISHED: December 21, 2017)

        raw: Vec<u8>
    }

    #[allow(dead_code)]
    pub enum Mode {
        CatA,
        CatB,
    }

    #[allow(dead_code)]
    pub enum Direction {
        Uplink,
        Downlink,
        Any,
    }

    #[allow(dead_code)]
    pub enum Suffix {
        Multiblock,
        Singleblock,
    }

    #[allow(dead_code)]
    impl Block {

        pub fn new(buffer: &[u8]) -> Result<Block, String>{

            Block::check_block(Block::from_slice(buffer))

        }

        fn from_slice(buffer: &[u8]) -> Block {

            let mut raw = Vec::new();
            
            for &byte in buffer.iter() {
                raw.push(match byte {
                    0x80        => NUL,
                    0x00..=0x80 => byte,
                    any_other   => Block::strip_parity_bit(any_other)
                })
            }

            Block { raw: raw }
        }

        pub fn from_unicode(unicode_string: String) -> Result<Block, String> {

            let mut raw = Vec::new();

            for byte in unicode_string.chars() {
                raw.push(match byte as u8 {
                    0x80 => NUL,
                    0x00..=0x80 => byte as u8,
                    any_other   => Block::strip_parity_bit(any_other)
                })
            }

            Block::check_block(Block { raw: raw })
        }

        fn check_block(acars_block: Block) -> Result<Block, String> {

            acars_block.get_mode() ?;

            acars_block.get_aircraft_array() ?;

            acars_block.get_ack() ?;

            acars_block.get_label_array() ?;

            acars_block.get_blk() ?;

            acars_block.get_stx() ?;

            acars_block.get_suffix() ?;

            Ok(acars_block)
        }

        fn strip_parity_bit(byte: u8) -> u8 {

            byte & 0b_0111_1111
        
        }

        pub fn clone(&self) -> Block {

            let clone  = self.raw.clone();

            Block { raw: clone }
        }

        pub fn get_mode(&self) -> Result<(u8, Mode), String> {

            let mode = match self.raw.get(0) {
                Some(mode) => mode,
                None =>  return Err(format!("Error getting Mode Character from Block [{:?}]", self.raw)),
            };

            match *mode {
                x if x == '2' as u8                   => Ok((x, Mode::CatA)),  // Uplink & Downlink (Cat A)
                x if x >= '`' as u8 && x <= '}' as u8 => Ok((x, Mode::CatB)),  // Uplink (Cat B)
                x if x >= '@' as u8 && x <= ']' as u8 => Ok((x, Mode::CatB)),  // Downlink (Cat B)
                x => { println!("{:?}", &self.raw); ; return Err(format!("Invalid Mode Character: {:#04x}", x)) } ,
            }
        }

        pub fn get_aircraft(&self) -> String {

            let mut aircraft = String::from("");

            match self.get_aircraft_array() { 
                
                Ok(result) => {
                    for &character in result.iter() {
                        match character {
                            x if x == NUL        => aircraft.push_str(""),
                            x if x == '.' as u8  => aircraft.push_str(""),
                            x                    => aircraft.push(x as char),
                        }
                    }
                },
                Err(_) => (),
            }

            aircraft
        }

        pub fn get_aircraft_array(&self) -> Result<[u8; 7], String> {

            let mut result = [NUL; 7];

            let aircraft_array = match self.raw.get(1..8) {
                Some(aircraft) => aircraft,
                None =>  return Err(format!("Error getting aircraft array from Block [{:?}]", self.raw)),
            };

            for (index, &value) in aircraft_array.iter().enumerate() {
                result[index] = match value {
                    x if x == NUL             => x,
                    x if x >= 0x20 && x < DEL => x,
                    x => return Err(format!("Invalid Aircraft Character: {:#04x}", x)),
                };
            }

            Ok(result)
        }

        pub fn get_ack(&self) -> Result<(u8, Direction), String> {

            let ack = match self.raw.get(8) {
                Some(ack) => ack,
                None =>  return Err(format!("Error getting Technical Acknowledgement Character from Block [{:?}]", self.raw)),
            };

            match *ack {
                x if x == NAK                         => Ok((x, Direction::Any)),
                x if x >= 'A' as u8 && x <= 'Z' as u8 => Ok((x, Direction::Downlink)),
                x if x >= 'a' as u8 && x <= 'a' as u8 => Ok((x, Direction::Downlink)),
                x if x >= '0' as u8 && x <= '9' as u8 => Ok((x, Direction::Uplink)),
                x => return Err(format!("Invalid Technical Acknowledgement Character: {:#04x}", x)),
            }
        }

        pub fn get_label(&self) -> String {
                
            let mut label = String::from("");

            match self.get_label_array() { 
                
                Ok(result) => {
                    for &character in result.iter() {
                        match character {
                            x if x == DEL => label.push_str("d"),
                            x             => label.push(x as char),
                        }
                    }
                },
                Err(_) => (),
            }

            label
        }

        pub fn get_label_array(&self) -> Result<[u8; 2], String> {

            let mut result = [NUL; 2];

            let label_array = match self.raw.get(9..11) {
                Some(label) => label,
                None =>  return Err(format!("Error getting Label array from Block [{:?}]", self.raw)),
            };

            for (index, &value) in label_array.iter().enumerate() {
                result[index] = match value {
                    x if x == DEL  => 'd' as u8,
                    x if x >= 0x20 => x,
                    x => return Err(format!("Invalid Label Character: {:#04x}", x)),
                };
            }

            Ok(result)
        }

        pub fn get_blk(&self) -> Result<(u8, Direction), String> {

            let blk = match self.raw.get(11) {
                Some(blk) => blk,
                None =>  return Err(format!("Error getting Block Identifier Character from Block [{:?}]", self.raw)),
            };

            match *blk { 
                x if x == NUL                         => Ok((x, Direction::Uplink)),
                x if x >= 'A' as u8 && x <= 'Z' as u8 => Ok((x, Direction::Uplink)),
                x if x >= 'a' as u8 && x <= 'a' as u8 => Ok((x, Direction::Uplink)),
                x if x >= '0' as u8 && x <= '9' as u8 => Ok((x, Direction::Downlink)),
                x => Err(format!("Invalid Block Identifier Character: {:#04x}", x)),
            }
        }

        pub fn get_stx(&self) -> Result<bool, String> {

            match self.raw.get(12) {
                Some(&STX) => Ok(true),
                Some(_)    => Ok(false),
                None       => Err(format!("Could not retrieve character at index 12")),
            }
        }

        pub fn get_suffix(&self) -> Result<(u8, usize, Suffix), String> {

            match self.raw.iter().position(|&x| (x == ETX || x == ETB)) {
                Some(index) => match index {
                    x if x >= 12 => match self.raw[x] {
                        ETB => Ok((self.raw[x], x, Suffix::Multiblock)),
                        ETX => Ok((self.raw[x], x, Suffix::Singleblock)),
                        any_other       => Err(format!("Invalid Suffix Character: {:#04x}", any_other)),
                    }
                    _                   => Err(format!("Block too small")),
                }
                None                    => Err(format!("Suffix not found")),
            }
        }

        pub fn get_raw(&self) -> String {

            let (_, size, _) = self.get_suffix().expect("Error at pub fn get_raw(): let (_, size, _) = self.get_suffix()");
            let mut raw = String::from("");

            for (index, &byte) in self.raw.get(..size)
                                          .expect("Error at pub fn get_raw(): for (index, &byte) in self.raw.get(..size)")
                                          .iter().enumerate() {
                let character = match byte {
                    x if x == DEL && index == 10 => 'd',
                    x if x == LF || x == CR      =>  x as char,
                    x if x <= 0x20 => { if index > 13                 { ' ' }
                        else                                          { '.' }
                    }, x                                             =>  x as char,
                };
                raw.push(character);
            }
            raw
        }

        pub fn get_essential(&self) -> String {

            let (_, direction) = self.get_blk().expect("Error at pub fn get_essential(): let (_, direction) = self.get_blk()");
            let (_, _, suffix) = self.get_suffix().expect("Error at pub fn get_essential(): let (_, _, suffix) = self.get_suffix()");

            let text = format!(" {} {}", 
                match direction {
                    Direction::Downlink => "↘",
                    Direction::Uplink   => "↗",
                    Direction::Any      => "↕",
                },
                match suffix {
                    Suffix::Multiblock  => "⊕",
                    Suffix::Singleblock => "⊝",
                },
            );

            let raw = self.get_raw();

            format!("{} [{}]", text, raw)       
        }
    }

    impl fmt::Display for Block {

        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {

            write!(f, "Block: [{}]", self.get_raw())
        }
    }
}

pub mod poa {
    
    extern crate chrono;

    use std::f64::consts::PI;
    use std::sync::Arc;
    use std::sync::mpsc;
    use chrono::{DateTime, Utc};
    use num::Complex;

    use super::ascii::*;
    use super::common::{Block, Reception};

    pub const CHANNEL_RATE: u32 = 12_500;


    #[allow(dead_code)]
    pub struct Channel {
        index: usize,
        central_frequency: f64,
        channel_frequency: f64,
        sample_rate: u32,
        sample_size: usize,
        channel_rate: u32,
        decimation_factor: u32,
        local_oscilator: Vec<Complex<f64>>,
        input: mpsc::Receiver<Arc<Vec<Complex<f64>>>>,
        msk: Demod,
        frame: Frame,
        output: mpsc::Sender<Reception>,
    }

    impl Channel {

        pub fn new(channel_steup: (usize, f64, f64, u32, usize, mpsc::Receiver<Arc<Vec<Complex<f64>>>>, mpsc::Sender<Reception>)) {

            let channel = Channel::setup(channel_steup);

            channel.process();
        }

        fn setup(channel_steup: (usize, f64, f64, u32, usize, mpsc::Receiver<Arc<Vec<Complex<f64>>>>, mpsc::Sender<Reception>)) 
            -> Channel {

            let (index, channel_frequency, central_frequency, sample_rate, sample_size, input, output) = channel_steup;

            let channel_rate = CHANNEL_RATE;
            let decimation_factor = sample_rate / channel_rate;

            let local_oscilator = Channel::local_oscilator(central_frequency, channel_frequency, sample_rate, sample_size);

            let msk = Demod::new(channel_rate);

            let frame = Frame::new();
 
            let channel = Channel {
                index,
                channel_frequency,
                central_frequency,
                sample_rate,
                channel_rate,
                decimation_factor,
                sample_size,
                local_oscilator,
                input,
                msk,
                frame,
                output,
            };

            eprintln!(" - ACARS decoder channel: {:.3} MHz", channel.channel_frequency);

            channel
        }

        fn process(mut self) {

            loop {

                let input = self.input.recv().expect("Error channel.input.recv()");

                let mut mixed_signal = Vec::with_capacity(input.len());
                let mut real_signal  = Vec::with_capacity(input.len() / self.decimation_factor as usize + 1);
                let mut local = Complex::new(0.0, 0.0);

                for (i, complex_value) in input.iter().enumerate() {


                    let mixed = complex_value * self.local_oscilator[i];

                    mixed_signal.push(mixed);

                    if i as u32 % self.decimation_factor == 0 {

                        let complex = local / self.decimation_factor as f64;
                        real_signal.push(complex.norm());
                        local = Complex::new(0.0, 0.0);

                    } else {

                        local += mixed;
                    }
                }

                self.demodulate(real_signal);
            }
        }

        fn demodulate(&mut self, real_signal: Vec<f64>) {

            let flen = (self.channel_rate / 1_200) + 1;
            let pllc = 3.8e-3;

            let mut idx = self.msk.idx;
            let mut o: i32;

            for i in 0..real_signal.len() {

                let s = 1800.0 / self.channel_rate as f64 * 2.0 * PI + self.msk.df;
                self.msk.phi += s;
                if self.msk.phi >= 2.0 * PI  { self.msk.phi -= 2.0 * PI; }

                let input = real_signal[i];
                self.msk.inb[idx as usize] = Complex::new(
                    input *  self.msk.phi.cos(),
                    input * -self.msk.phi.sin()
                );
                idx = (idx + 1) % flen as i32;

                self.msk.clk += s;
                if self.msk.clk >= 3.0 * PI / 2.0 {
                    self.msk.clk -= 3.0 * PI / 2.0;

                    o = flen as i32 - idx;
                    let mut v = Complex::new(0.0, 0.0);
                    for j in 0..flen {
                        v += self.msk.h[o as usize] * self.msk.inb[j as usize];
                        o += 1;
                    }

                    let lvl = v.norm();
                    v /= lvl + 1e-6;
                    self.msk.lvl = 0.99 * self.msk.lvl + 0.01 * lvl/5.2;

                    let (dphi, bit) = Demod::quadrature(self.msk.s, v);

                    self.decode(bit);

                    self.msk.s += 1;
                    self.msk.df = pllc * dphi;
                }
            }

            self.msk.idx = idx;
        }

        fn decode(&mut self, bit: bool) {

            self.frame.byte >>= 1;
            if bit == true { self.frame.byte |= 0b_1000_0000; }
            let no_parity =  self.frame.byte  & 0b_0111_1111;
            self.frame.remaining_bits -= 1;

            if self.frame.remaining_bits <= 0 {
                match self.frame.stage {

                    FrameStage::WSYN1 => {

                        if no_parity == SYN { 
                            self.frame.stage = FrameStage::WSYN2; 
                            self.frame.remaining_bits = 8;
                        } else if no_parity == !SYN {
                            self.frame.stage = FrameStage::WSYN2; 
                            self.frame.remaining_bits = 8;
                            self.msk.s ^= 0b_0010;
                        } 
                    },

                    FrameStage::WSYN2 => {

                        if no_parity == SYN {
                            self.frame.stage = FrameStage::WSOH;
                            self.frame.remaining_bits = 8; 
                        } else if no_parity == !SYN {
                            self.frame.stage = FrameStage::WSYN2; 
                            self.frame.remaining_bits = 8;
                            self.msk.s ^= 0b_0010;
                        } else {
                            self.frame.clear(); 
                            self.msk.df = 0.0;
                        }
                    },

                    FrameStage::WSOH => {

                        if no_parity == SOH {

                            // println!("SOH");
                            self.frame.stage = FrameStage::BLK;
                            self.frame.remaining_bits = 8;

                            self.frame.start = Utc::now();
                            self.msk.lvl = 0.0;

                        } else {
                            // println!("--H");
                            self.frame.clear();
                        }
                    },

                    FrameStage::BLK => {
                        
                        self.frame.bytes.push(self.frame.byte);

                        if ! Frame::is_parity_correct(&self.frame.byte) {
                            self.frame.parity_errors += 1;
                        }

                        if self.frame.parity_errors >= 13 {

                            println!("Too much parity errors");
                            // println!("{:?}", &self.frame.bytes);
                            for each_byte in &self.frame.bytes {
                                print!("{}", (*each_byte & 0b_0111_1111) as char);
                            }
                            println!(" ");

                            self.frame.clear();
                        }

                        if no_parity == ETX || no_parity == ETB { 

                            self.frame.stage = FrameStage::CRC1; 
                            self.frame.remaining_bits = 8;
                            self.frame.end = Utc::now();

                            if self.frame.bytes.len() < 13 {

                                println!("Too small");
                                self.frame.clear();
                            } 

                        } else { 
                            self.frame.remaining_bits = 8;
                        }
                            
                    },

                    FrameStage::CRC1 => {

                        // println!("CRC1");
                        self.frame.crc[0] = self.frame.byte;
                        self.frame.stage = FrameStage::CRC2;
                        self.frame.remaining_bits = 8;
                        
                    },

                    FrameStage::CRC2 => {

                        // println!("CRC2");
                        self.frame.crc[1] = self.frame.byte;
                        self.frame.stage = FrameStage::END;
                        self.frame.remaining_bits = 8;
                    },

                    FrameStage::END => {

                        Channel::build_output(&self);

                        self.frame.clear();

                        if  no_parity != DEL {
                            // println!("NO END");
                        }

                        self.msk.df = 0.0;
                    },

                }
            }
        }

        fn build_output(channel: &Channel) {

            // ERROR CHECKING & CORRECTION
            // self.frame.check_and_correct();

            // CHECK CRC

            // PARITY CHECK

            let level = (10.0 * channel.msk.lvl.log10()) as isize;
            let station = String::from("TEST");

            let no_parity = Frame::remove_parity(&channel.frame.bytes);
            let block = match Block::new(no_parity.as_slice()) {
                Ok(block) => Some(block),
                Err(error) => {
                    println!("ERROR building block: {}\n{:?}", error, no_parity);
                    None
                }
            };

            let reception = Reception {
                date_time: channel.frame.start,
                station: station,
                channel: channel.index,
                frequency: channel.channel_frequency,
                level: level,
                acars_block: block
            };     

            channel.output.send(reception).expect("failure sending block");
        }

        fn local_oscilator(central_frequency: f64, channel_frequency: f64, sample_rate: u32, sample_size: usize) -> Vec<num::Complex<f64>> {

            let freq_difference = channel_frequency - central_frequency;

            let shift_factor = freq_difference * 1e6 / sample_rate as f64;

            let mut local_oscilator: Vec<num::Complex<f64>> = Vec::new();

            for i in 0..sample_size * 1024 {
                let shift = -shift_factor * 2.0 * PI * i as f64;
                local_oscilator.push(Complex::from_polar(&1.0, &shift));
            }

            local_oscilator
        }
    }

    #[derive(Debug)]
    struct Demod {
        channel_rate: u32,
        df:  f64,
        phi: f64,
        clk: f64,
        lvl: f64,
        s:   u32,
        idx: i32,
        h:   Vec<f64>,
        inb: Vec<Complex<f64>>,
    }

    impl Demod {

        pub fn new(channel_rate: u32) -> Demod {

            let flen = (channel_rate as usize / 1_200) + 1;

            let mut h = vec![0.0; 2*flen];
            let inb = vec![Complex::new(0.0, 0.0); flen];

            for i in 0..flen {

                let value = (2.0 * PI * 600.0 / channel_rate as f64 * (i as i32 - flen as i32/2) as f64).cos();

                h[i]      = value;
                h[i+flen] = value;
            }
           
            Demod {
                channel_rate,
                df:  0.0,
                phi: 0.0,
                clk: 0.0,
                lvl: 0.0,
                s:   0,
                idx: 0,
                h:   h,
                inb: inb,
            }
        }

        fn quadrature(msk_s: u32, v: Complex<f64>) -> (f64, bool) {

            let dphi: f64;

            let vo = match msk_s & 3 {
                0 => { if v.re >= 0.0 { dphi =  v.im; } else { dphi = -v.im; } ;  v.re },
                1 => { if v.im >= 0.0 { dphi = -v.re; } else { dphi =  v.re; } ;  v.im },
                2 => { if v.re >= 0.0 { dphi =  v.im; } else { dphi = -v.im; } ; -v.re },
                3 => { if v.im >= 0.0 { dphi = -v.re; } else { dphi =  v.re; } ; -v.im },
                _ => panic!("Invalid result of match self.s & 3"),
            };

            let bit = match vo {
                x if x >  0.0 => true,
                x if x <= 0.0 => false,
                _             => panic!("Invalid vo"),
            };

            (dphi, bit)
        }
    }

    #[allow(dead_code)]
    #[derive(Debug)]
    struct Frame {
        byte: u8,
        stage: FrameStage,
        remaining_bits: i32,
        bytes: Vec<u8>,
        size: usize,
        start: DateTime<Utc>,
        end: DateTime<Utc>,
        lvl_vec: Vec<f64>,
        parity_errors: usize,
        crc: [u8; 2],
    }

    #[allow(dead_code)]
    #[derive(Debug)]
    enum FrameStage {
        WSYN1,
        WSYN2,
        WSOH,
        BLK,
        CRC1,
        CRC2,
        END,
    }

    impl Frame {

        pub fn new() -> Frame {

            let stage = FrameStage::WSYN1;
            let byte = 0b_0000_0000;
            let remaining_bits = 1;
            let bytes = Vec::new();
            let size = 0;
            let parity_errors = 0;
            let crc  = [0; 2];

            let now: chrono::DateTime<Utc> = Utc::now();

            Frame {
                byte,
                stage,
                remaining_bits,
                bytes,
                size,
                start: now,
                end: now,
                lvl_vec: Vec::new(),
                parity_errors,
                crc,
            }
        }

        pub fn is_parity_correct(byte: &u8) -> bool {

            let mut bit_count = 0;

            for i in 0..8 {
                    
                let mask = 0b_0000_0001 * 2_u8.pow(i as u32);
                        
                if byte & mask == mask {
                            
                    bit_count += 1;
                            
                }
                        
                // if print == true { println!("{:0>8b} & {:0>8b} = {:0>8b} ({}) ({})", byte, mask, byte & mask, bit_count, bit_count % 2 == 1); }
            }

            bit_count % 2 == 1

        }

        pub fn remove_parity(bytes: &Vec<u8>) -> Vec<u8> {

            let mut no_parity: Vec<u8> = Vec::new();

            for byte in bytes {

                no_parity.push(byte & 0b_0111_1111 );

            }

            no_parity
        }

        pub fn clear(&mut self) {

            self.stage = FrameStage::WSYN1;
            self.byte = 0b_0000_0000;
            self.remaining_bits = 1;
            self.bytes = Vec::new();
            self.size = 0;
            self.parity_errors = 0;
            self.crc  = [0; 2];

        }
    }
}

pub mod vdl {

    // to be defined

}

#[allow(dead_code)]
pub mod ascii {

    pub const NUL: u8 = 0x00;
    pub const SOH: u8 = 0x01;
    pub const STX: u8 = 0x02;
    pub const ETX: u8 = 0x03;
    pub const EOT: u8 = 0x04;
    pub const ENQ: u8 = 0x05;
    pub const ACK: u8 = 0x06;
    pub const BEL: u8 = 0x07;
    pub const BS:  u8 = 0x08;
    pub const HT:  u8 = 0x09;
    pub const LF:  u8 = 0x0a;
    pub const VT:  u8 = 0x0b;
    pub const FF:  u8 = 0x0c;
    pub const CR:  u8 = 0x0d;
    pub const SO:  u8 = 0x0e;
    pub const SI:  u8 = 0x0f;

    pub const DLE: u8 = 0x10;
    pub const DC1: u8 = 0x11;
    pub const DC2: u8 = 0x12;
    pub const DC3: u8 = 0x13;
    pub const DC4: u8 = 0x14;
    pub const NAK: u8 = 0x15;
    pub const SYN: u8 = 0x16;
    pub const ETB: u8 = 0x17;
    pub const CAN: u8 = 0x18;
    pub const EM:  u8 = 0x19;
    pub const SUB: u8 = 0x1a;
    pub const ESC: u8 = 0x1b;
    pub const FC:  u8 = 0x1c;
    pub const GS:  u8 = 0x1d;
    pub const RS:  u8 = 0x1e;
    pub const DEL: u8 = 0x7f;
}

pub mod crc {   

    // to be defined

}
