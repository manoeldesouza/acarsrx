/*
 *
 *  This code is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License version 2
 *  published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

//! # acarsrx
//! Copyright (c) 2019 Manoel Souza <manoel.desouza@outlook.com.br>
//! 
//! ACARS decoder for RTL-SDR and (in future) other SDR devices
//! 
//! Inspired on Acarsdec Copyright (c) by Thierry Thierry Leconte
//! 
//! https://github.com/TLeconte/acarsdec
//!
//! and Acars-sdrplay Copyright (c) by Jan van Katwijk <J.vanKatwijk@gmail.com>
//! 
//! https://github.com/JvanKatwijk/acars-sdrplay
//! 
//! 
//!  ## Introduction:
//!
//! acarsrx is designed to be able to run each SDR device in separate thread enabling multiple devices 
//! to be controlled from a single execution. For each SDR device controlled, multiple channels can be 
//! used to demodulate and decode multiple frequencies. Special attention is put into decoupling the 
//! SDR handling and the channels processing (where demodulation, decoding and block formating). This 
//! way, it is planned for other SDR devices to be controlled from the same application (HackRF is a 
//! work in progress). Finally, acarsrx is designed in a way to recover from SDR devices disconnects.
//! Being USB devices the author experienced common disconnect situations and enable the application 
//! to recover from such situation proved to be an asset. 
//! 
//! 
//! ## Design:
//! 
//! Each SDR device is planned to have all instructions handling from a single Rust file (as done for 
//! rtlsdr.rs and being done for hackrf.rs).
//! 
//! Following the same concept, each protocol is planned to be fully defined in a single Rust file (as 
//! done for acars.rs). acars.rs defines all instructions for plain old ACARS (POA) and is planned to 
//! also include VDL Mode 2 decoding in future releases.
//!
//! The system is designed to initiate threads to handle each SDR device and other threads for each 
//! channel from within the SDR thread. 
//! 
//! 
//! ### Thread execution example: 
//! 
//! ```bash
//! acarsrx --rtlsdr 0 131.550 131.725 --rtlsdr 1 129.125 129.425
//! ```
//! 
//! ```text
//!                         / acars.rs 131.550 \
//!           -- rtlsdr.rs 0                    \
//!         /               \ acars.rs 131.725 \ \
//! main.rs                                      output.rs
//!         \               / acars.rs 129.125 / /
//!           -- rtlsdr.rs 1                    /
//!                         \ acars.rs 129.425 /
//! ```
//!
//! ### Sample Output:
//! ```text
//! 07:41:55 RADIO-r0-0 131.550 MHz -24 ↗ ⊝ [2........SQ..01XAYULCYUL1ARINC]
//! 07:42:07 RADIO-r0-0 129.125 MHz -21 ↗ ⊝ [2.C-FTJQ1_dY]
//! 07:42:39 RADIO-r0-0 131.725 MHz -14 ↘ ⊝ [2..N8767K_d5.S75AZD8767]
//! 07:45:11 RADIO-r0-0 131.725 MHz -23 ↗ ⊝ [2........SQ..02XSYULCYUL04527N07344WV136975/]
//! 07:46:05 RADIO-r0-0 131.550 MHz -21 ↗ ⊝ [2........SQ..02XAYULCYUL14528N07344WV136975/ARINC]
//! 07:46:36 RADIO-r0-0 131.725 MHz -20 ↘ ⊝ [2..N8767N_d1.S78AZD8767]
//! -------- ---------- ----------- --- - - -------------------------
//!     |       |          |         |  | | ACARS Block
//!     |       |          |         |  | SingleBlock or MultiBlock
//!     |       |          |         |  Uplink or Downlink
//!     |       |          |         Signal Strenght
//!     |       |          Channel Frequency
//!     |       Station, device & Channel number 
//!     Time (hh:mm:ss)
//! ``` 
//! 
//! 
//! ## Challenges:
//! 
//! There is some odd behavior still under investigation by the author. The simple sample_size achived by 
//! diving the device SAMPLE_RATE per the CHANNEL_RATE, shows poor decoding performance, with improvements 
//! being a achived with multiple of these values. To achieve better performance the author applied a 
//! factor constant named SIZE_RATE_MULTIPLIER = 20. Values from 10 in general are enough for good results.
//! When compared with other acars decoders this seems to be a deffect still to be investigated.
//! 
//! 
//! ## Next steps:
//! 
//! There are some functionalities still not present in acarsrx to achieve V1.0, like:
//!
//! 1. Multiple Reception qualities to be submitted to output.rs (Bad CRC block, Bad parity bit, etc). It 
//!    will be for the output.rs thread to define what is to be printed and how.
//! 
//! 2. CRC error checking, correction and re parity check.
//! 
//! Other improvements, new functionalities are planned for future:
//! 
//! 3. Implement multiple output formats in addition to current single-line print, for example JSON, Ncurses
//!    based interface, UDP transmission, Save to file. 
//! 
//! 4. Enable ACARS Block assembly into ACARS Messages (Type-B messages) with translation of Labels into 
//!    SMI, and standard message formating (Following AEEC 620-8 specification).
//! 
//! 5. Control of HackRF SDR. Other SDR devices can be controled as well given that similar methods to
//!    enable the IQ buffer to be retrieved and transformed into vectors of complex numbers.
//! 
//! 6. Implementation of VDL Mode 2 decoding. By design a single SDR should enable decoding of multiple 
//!    protocols and both Plain Old ACARS and VDL Mode 2 should be possible to be decoded, give that enough
//!    bandwidth and channel separation is provided.
//! 
//! 
//! ## Future developments:
//! 
//! The exercise on implementing a decoupling between ACARS demodulation and decoding from the sample
//! extraction from SDR device, inspires the author to consider developing a general SDR driver which
//! would be able to control and generate samples for multiple protocols (couting on community support)
//! and outputs. Such server would enable multiple decoders to consume the same sources and multiple 
//! outputs to be triggered from those.
//! 
//! | Input |  Process  |  Output |
//! |:-----:|:---------:|:-------:|
//! |RTLSDR | AM FM SSB | speaker |
//! |HACKRF | ACARS     | wav     |
//! |AIRSPY | CW        | json    |
//! |etc    | APRS      | txt     |
//! |       | etc       | tcp/udp |
//! 
//! Suggested name: sdrD (software defined radio Daemon).
//! 
//! ## Additional details on acarsrx:
//! 
//! The modules section below describes the internal functions in acarsrx
//! 
//! My regards, Manoel Souza
//! 08-Oct-2019
//! 
use std::env;
use std::fs;
use std::thread;
use std::sync::mpsc;

#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate serde_json;

extern crate chrono;

mod rtlsdr;
mod hackrf;
mod acars;
mod output;


fn main() {

    banner();

    let args: Vec<String> = env::args().collect();
    let config = Config::from_arguments(args);
    let radios = config.radios.len();
    let (tx, rx) = mpsc::channel();

    for radio in config.radios {

        let number = radio.number;
        let model = radio.model;
        let frequencies = radio.frequencies.clone();
        let tx_local = tx.clone();

        match model {
            Model::RTLSDR => { thread::spawn(move || { rtlsdr::Device::new(number, frequencies, tx_local); }); },
            Model::HACKRF => { thread::spawn(move || { hackrf::Device::new(number, frequencies, tx_local); }); },
        }
    }
    
    if radios > 0 { 

        let output = thread::spawn(|| { output::thread(rx); } );
        output.join().unwrap(); 
    
    } else { 

        usage();
    }
}

fn banner() {

    let version = String::from("v0.8.0");
    let release = String::from("07-Oct-2019");

    println!("
 acarsrx {} - {} 
 Copyright (c) 2019 Manoel Souza <manoel.desouza@outlook.com.br>
    ", version, release);
}

fn usage() {

    let default_config = String::from(
            "{\n\t\"forward_to\": \"\",\n\t\"save_to\": \"\",\n\t\"radios\": [\n\
            \t\t{ \"model\": \"RTLSDR\", \"number\": 0, \"mode\": \"ACARS\", \
            \"frequencies\": [ 131.550, 131.725 ] }\n\t]\n}\n"
    );

    println!(" Usage: acarsrx [ --rtlsdr index freq1 freq2 ... freqN ] [ --config filename ]
        \
        \n Ex: acarsrx --rtlsdr 0 131.550 131.725 --rtlsdr 1 129.125
        \
        \n --rtlsdr: Controls a single RTL-SDR dongle.\
        \n   * index: RTL-SDR device id found with rtl_test program right after \"Found X device(s):\"\
        \n     Multiple rtl-sdr devices can be operated in a single execution.\
        \n   * freqN: List of frequencies to be decoded. Make sure each RTL-SDR device is set with \
        \n     frequencies no farther than 1 MHz.
        \
        \n --config: Instead of entering all parameters via command line, JSON formated configuration \
        \n           file is used instead. Default contents are:\
        \n
        \nDefault Config:\
        \n{}\
        \nSample Output:
        \
        \n 07:41:55 RADIO-r0-0 131.550 MHz -24 ↗ ⊝ [2........SQ..01XAYULCYUL1ARINC]
        \n 07:42:07 RADIO-r0-0 129.125 MHz -21 ↗ ⊝ [2.C-FTJQ1_dY]
        \n 07:42:39 RADIO-r0-0 131.725 MHz -14 ↘ ⊝ [2..N8767K_d5.S75AZD8767]
        \n 07:45:11 RADIO-r0-0 131.725 MHz -23 ↗ ⊝ [2........SQ..02XSYULCYUL04527N07344WV136975/]
        \n 07:46:05 RADIO-r0-0 131.550 MHz -21 ↗ ⊝ [2........SQ..02XAYULCYUL14528N07344WV136975/ARINC]
        \n 07:46:36 RADIO-r0-0 131.725 MHz -20 ↘ ⊝ [2..N8767N_d1.S78AZD8767]
        \n -------- ---------- ----------- --- - - -------------------------
        \n     |       |          |         |  | | ACARS Block
        \n     |       |          |         |  | SingleBlock or MultiBlock
        \n     |       |          |         |  Uplink or Downlink
        \n     |       |          |         Signal Strenght
        \n     |       |          Channel Frequency
        \n     |       Station, device & Channel number 
        \n     Time (hh:mm:ss)
        \n", default_config);

}

#[derive(Serialize, Deserialize, Debug)]
pub struct Config {
    pub radios: Vec<Radio>,
    pub forward_to: String,
    pub save_to: String,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Radio {
    pub model: Model,
    pub number: i32,
    pub mode: Mode,
    pub frequencies: Vec<f64>,
}

#[derive(Serialize, Deserialize, Debug)]
pub enum Model {
    RTLSDR,
    HACKRF,
}

#[derive(Serialize, Deserialize, Debug)]
pub enum Mode {
    ACARS,
    VDLM2,
}

impl Config {

    pub fn load(config_file: &String) -> Config {

        let default_config = String::from(
            "{\n\t\"forward_to\": \"\",\n\t\"save_to\": \"\",\n\t\"radios\": [\n\
            \t\t{ \"model\": \"RTLSDR\", \"number\": 0, \"mode\": \"ACARS\", \"frequencies\": [ 131.550, 131.725 ] }\n\t]\n}\n"
        );

        let serialized = match fs::read_to_string(config_file) {
            Ok(file_contents) => file_contents,
            Err(_) => {
                eprintln!("Error opening config.json. Applying default config");
                Config::save(&default_config, config_file);
                default_config 
            },
        };    

        let deserialized: Config = serde_json::from_str(&serialized).expect(
            "Error reading config file. Delete config.json and a default configuration will be generated."
        );

        deserialized
    }

    pub fn save(config: &String, config_file: &String) {

        match fs::write(config_file, config) {
            Ok(_) => (),
            Err(_) => eprintln!("Error saving default config"),
        }

    }

    pub fn from_arguments(args: Vec<String>) -> Config {

        let mut forward_to =  String::from("");
        let mut save_to =  String::from("");
        let mut radios = Vec::new();

        let mut i  = 1;
        while i < args.len() {

            if &args[i] == "--config" {

                let config_file = String::from(args[i+1].trim());
                let config = Config::load(&config_file);

                return config;

            } else if &args[i] == "--save_to" {

                save_to.push_str(args[i+1].trim());
                i += 2;

            } else if &args[i] == "--forward_to" {

                forward_to.push_str(args[i+1].trim());
                i += 2;
            
            } else if &args[i] == "--rtlsdr" {

                let model = Model::RTLSDR;
                let mode = Mode::ACARS;
                let number = args[i+1].trim().parse().unwrap();

                radios.push(Radio { 
                    model,
                    mode,
                    number,
                    frequencies: Vec::new(),
                });

                i += 2;

            } else if &args[i] == "--hackrf" {

                let model = Model::HACKRF;
                let mode = Mode::ACARS;
                let number = args[i+1].trim().parse().unwrap();

                radios.push(Radio { 
                    model,
                    mode,
                    number,
                    frequencies: Vec::new(),
                });

                i += 2;
            
            } else {

                let frequency: f64 = args[i].trim().parse().unwrap();
                let this_radio = &radios.len()-1;

                radios[this_radio].frequencies.push(frequency);
                i += 1;
            }
        }

        Config {
            save_to,
            forward_to,
            radios,
        }
    }
}
