
use std::f64::consts::PI;
use std::str;
use std::sync::Arc;
use std::sync::mpsc;
use num::Complex;


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
}

impl Channel {

    pub fn new(channel_steup: (usize, f64, f64, u32, usize, mpsc::Receiver<Arc<Vec<Complex<f64>>>>)) {

        let channel = Channel::setup(channel_steup);

        Channel::process(channel);
    }

    fn setup(channel_steup: (usize, f64, f64, u32, usize, mpsc::Receiver<Arc<Vec<Complex<f64>>>>)) -> Channel {

        let (index, channel_frequency, central_frequency, sample_rate, sample_size, input) = channel_steup;

        let channel_rate = CHANNEL_RATE;
        let decimation_factor = sample_rate / channel_rate;

        let local_oscilator = Channel::local_oscilator(central_frequency, channel_frequency, sample_rate, sample_size);

        let msk = Demod::new(channel_rate);

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
            msk
        };

        eprintln!(" - ACARS decoder channel: {:.3} MHz", channel.channel_frequency);

        channel
    }

    fn process(mut channel: Channel) {

        loop {

            let input = channel.input.recv().expect("Error channel.input.recv()");

            let mut mixed_signal = Vec::with_capacity(input.len());
            let mut real_signal  = Vec::with_capacity(input.len() / channel.decimation_factor as usize + 1);
            let mut local = Complex::new(0.0, 0.0);

            for (i, complex_value) in input.iter().enumerate() {

                let mixed = complex_value * channel.local_oscilator[i];

                mixed_signal.push(mixed);

                if i as u32 % channel.decimation_factor == 0 {

                    let complex = local / channel.decimation_factor as f64;
                    real_signal.push(complex.norm());
                    local = Complex::new(0.0, 0.0);

                } else {

                    local += mixed;
                }
            }

            channel.msk.demodulate(real_signal);
        }
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
    frame: Frame,
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

        let frame = Frame::new();
        
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
            frame: frame,
        }
    }

    pub fn demodulate(&mut self, real_signal: Vec<f64>) {

        let flen = (self.channel_rate / 1_200) + 1;
        let pllc = 3.8e-3;

        let mut idx = self.idx;
        let mut o: i32;

        for i in 0..real_signal.len() {

            let s = 1800.0 / self.channel_rate as f64 * 2.0 * PI + self.df;
            self.phi += s;
            if self.phi >= 2.0 * PI  { self.phi -= 2.0 * PI; }

            let input = real_signal[i];
            self.inb[idx as usize] = Complex::new(
                input *  self.phi.cos(),
                input * -self.phi.sin()
            );
            idx = (idx + 1) % flen as i32;

            self.clk += s;
            if self.clk >= 3.0 * PI / 2.0 {
                self.clk -= 3.0 * PI / 2.0;

                o = flen as i32 - idx;
                let mut v = Complex::new(0.0, 0.0);
                for j in 0..flen {
                    v += self.h[o as usize] * self.inb[j as usize];
                    o += 1;
                }

                let lvl = v.norm();
                v /= lvl + 1e-6;
                self.lvl = 0.99 * self.lvl + 0.01 * lvl/5.2;

                let (dphi, bit) = Demod::quadrature(self.s, v);

                self.frame.decode(bit, self.s, self.df);

                self.s += 1;
                self.df = pllc * dphi;
            }
        }

        self.idx = idx;
        // self.phi = p;
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
pub struct Frame {
    byte: u8,
    stage: FrameStage,
    remaining_bits: i32,
    acars_frame: Vec<u8>,
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
        let byte: u8 = 0b_0000_0000;
        let remaining_bits: i32 = 1;
        let acars_frame: Vec<u8> = Vec::new();

        Frame {
            byte,
            stage,
            remaining_bits,
            acars_frame,
        }
    }

    pub fn decode(&mut self, bit: bool, mut _msk_s: u32, mut _msk_df: f64) {

        const SOH: u8 = 0x01;
        const ETX: u8 = 0x03;
        const SYN: u8 = 0x16;
        const ETB: u8 = 0x17;
        const DEL: u8 = 0x7f;

        self.byte >>= 1;
        if bit == true { self.byte |= 0b_1000_0000; }
        let no_parity = self.byte & 0b_0111_1111;
        self.remaining_bits -= 1;

        if self.remaining_bits <= 0 {
            match self.stage {

                FrameStage::WSYN1 => {

                    if no_parity == SYN { 
                        self.stage = FrameStage::WSYN2; 
                        self.remaining_bits = 8;
                    } else if no_parity == !SYN {
                        self.stage = FrameStage::WSYN2; 
                        self.remaining_bits = 8;
                        _msk_s ^= 0b_0010;
                    } else {
                        self.remaining_bits = 1; }
                },

                FrameStage::WSYN2 => {

                    if no_parity == SYN {
                        println!("\nSYN2");
                        self.stage = FrameStage::WSOH;
                        self.remaining_bits = 8; 
                    } else if no_parity == !SYN {
                        self.stage = FrameStage::WSYN2; 
                        self.remaining_bits = 8;
                        _msk_s ^= 0b_0010;
                    } else { 
                        self.stage = FrameStage::WSYN1;
                        self.remaining_bits = 1; 
                        _msk_df = 0.0; }
                },

                FrameStage::WSOH => {

                    if no_parity == SOH {
                        println!("SOH");
                        self.stage = FrameStage::BLK;
                        self.remaining_bits = 8;
                        self.acars_frame.push(no_parity);
                    } else {
                        println!("--H");
                        self.stage = FrameStage::WSYN1;
                        self.remaining_bits = 1; }
                },

                FrameStage::BLK => {
                    
                    if  no_parity == ETX || no_parity == ETB { 
                        self.stage = FrameStage::CRC1; 
                        self.remaining_bits = 8;
                        println!("{}", str::from_utf8(&self.acars_frame).expect("FAILURE INTERPRETING BLOCK")); 
                        self.acars_frame.clear();
                    } else { 
                        self.remaining_bits = 8; }
                        self.acars_frame.push(no_parity);
                },

                FrameStage::CRC1 => {

                    println!("CRC1");
                    self.stage = FrameStage::CRC2;
                    self.remaining_bits = 8;
                },

                FrameStage::CRC2 => {

                    println!("CRC2");
                    self.stage = FrameStage::END;
                    self.remaining_bits = 8;
                },

                FrameStage::END => {

                    if  no_parity == DEL { 
                        println!("END");
                        self.stage = FrameStage::WSYN1;
                        self.remaining_bits = 1;
                    } else {
                        println!("NO END");
                        self.stage = FrameStage::WSYN1;
                        self.remaining_bits = 1;
                    }

                    _msk_df = 0.0;
                },

            }
        }
    }
}
