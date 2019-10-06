
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
mod acars;
mod output;

fn main() {

    let args: Vec<String> = env::args().collect();
    let config = Config::from_arguments(args);

    let (tx, rx) = mpsc::channel();
    let output = thread::spawn(|| { output::thread(rx); } );

    for radio in config.radios {

        let number = radio.number;
        let model = radio.model;
        let frequencies = radio.frequencies.clone();

        let tx_local = tx.clone();

        match model {
            Model::RTLSDR => { thread::spawn(move || { rtlsdr::Device::new(number, frequencies, tx_local); }); },
            Model::HACKRF => { },
        }
    }
    
    output.join().unwrap();
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
//    ADSB,
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
