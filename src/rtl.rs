
use std::thread;
use std::sync::Arc;
use std::sync::mpsc;
use num::Complex;
use crate::acars;

const SAMPLE_RATE: u32   = 2_400_000;


#[allow(dead_code)]
pub struct Device {
    dev: rtlsdr::RTLSDRDevice,
    index: i32,
    central_frequency: f64,
    sample_rate: u32,
    buffer_size: usize,
    outputs: Vec<mpsc::Sender<Arc<Vec<Complex<f64>>>>>,
}

impl Device {

    pub fn new(index: i32, frequencies: Vec<f64>) {

        loop {

            let device = match Device::setup(index, &frequencies) {
                Ok(device) => device,
                Err(_) => {
                    eprintln!("Restarting RTL-SDR device {} in 5 secs...", index);
                    std::thread::sleep(std::time::Duration::from_secs(5));
                    continue;
                },
            };

            Device::process(device);
        }
    }

    fn setup(index: i32, frequencies: &Vec<f64>) -> Result<Device, rtlsdr::RTLSDRError> {

        let central_frequency = Device::central_frequency(&frequencies);
        let sample_rate = SAMPLE_RATE;
        let sample_size = (sample_rate / acars::CHANNEL_RATE) as usize;
        let buffer_size = sample_size * 2 * 1024;

        let mut dev = rtlsdr::open(index) ?;

        dev.set_center_freq((central_frequency * 1e6) as u32) ?;
        dev.set_sample_rate(sample_rate) ?;
        dev.set_tuner_gain_mode(false) ?;
        dev.reset_buffer() ?;

        let sample_rate = dev.get_sample_rate() ?;
        let central_frequency = dev.get_center_freq() ? as f64 / 1e6;

        let mut device = Device {
            dev,
            index,
            central_frequency,
            sample_rate,
            buffer_size,
            outputs: Vec::new()
        };

        for (i, frequency) in frequencies.iter().enumerate() {
            let (tx, rx) = mpsc::channel();
            let channel_setup = (
                i,
                *frequency,
                central_frequency,
                sample_rate,
                sample_size,
                rx
            );

            thread::spawn(move || { acars::Channel::new(channel_setup); });

            device.outputs.push(tx);
        }


        eprintln!("RTL-SDR device {} central frequency: {:.3} MHz", index, central_frequency);
        Ok(device)
    }

    fn process(mut device: Device) {

        loop {

            let byte_buffer = match device.dev.read_sync(device.buffer_size) {
                Ok(buffer) => buffer,
                Err(_)     => break,
            };

            let complex_buffer = Arc::new({
                let mut complex_buffer: Vec<Complex<f64>> = Vec::with_capacity(device.buffer_size/2);
                for iq in byte_buffer.chunks(2) {

                    let r = iq[0] as f64 / (255.0/2.0) - 1.0;
                    let g = iq[1] as f64 / (255.0/2.0) - 1.0;
                    let complex = Complex::new(r, g);
                    
                    complex_buffer.push(complex);
                }
                complex_buffer
            });

            for output in device.outputs.iter() {
                output.send(complex_buffer.clone()).expect("Error output.send(complex_buffer.clone())");
            }

        }
    }

    fn central_frequency(frequencies: &Vec<f64>) -> f64{

        let mut center_frequency = 0.0;

        for frequency in frequencies.iter() {
            center_frequency += frequency;
        }

        center_frequency /= frequencies.len() as f64;

        if frequencies.len() == 1 {
            center_frequency += 0.025;
        }

        center_frequency
    }
}