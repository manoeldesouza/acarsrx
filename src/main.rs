
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

    let frequencies0 = vec!(
        // 129.125,
        131.550,
        131.725,
    );

    let frequencies1 = vec!(
        129.125,
    );

    let (tx, rx) = mpsc::channel();
    thread::spawn(|| { output::thread(rx); } );

    let tx0 = tx.clone();
    let tx1 = tx.clone();

    let handle = thread::spawn(|| { rtlsdr::Device::new(0, frequencies0, tx0); });
    let handle = thread::spawn(|| { rtlsdr::Device::new(1, frequencies1, tx1); });
    handle.join().unwrap();
}



#[derive(Serialize, Deserialize, Debug)]
pub struct Config {
    pub radios: Vec<Radio>,
    pub forward_to: String,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Radio {
    pub device_number: usize,
    pub mode: Mode,
    pub frequencies: Vec<f64>,
}

#[derive(Serialize, Deserialize, Debug)]
pub enum Mode {
    ACARS,
    VDLM2,
//    ADSB,
}

impl Config {

    pub fn load() -> Config {

        let config_file = "./config.json";
        let default_config = String::from(
            "{\n\t\"forward_to\": \"\",\n\t\"save_to\": \"\",\n\t\"radios\": [\n\
            \t\t{ \"device_number\": 0, \"mode\": \"ACARS\", \"frequencies\": [ 131.550, 131.725 ] }\n\t]\n}\n"
        );

        let serialized = match fs::read_to_string(config_file) {
            Ok(file_contents) => file_contents,
            Err(_) => {
                eprintln!("Error opening config.json. Applying default config");
                Config::save(&default_config);
                default_config 
            },
        };    

        let deserialized: Config = serde_json::from_str(&serialized).expect(
            "Error reading config file. Delete config.json and a default configuration will be generated."
        );

        deserialized
    }

    pub fn save(config: &String) {

        match fs::write("./config.json", config) {
            Ok(_) => (),
            Err(_) => eprintln!("Error saving default config"),
        }

    }
}
