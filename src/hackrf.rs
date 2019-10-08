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


//! hackrf.rs provides all functionalities for HACKRF SDR device handling. The Device object 
//! once instantiated will create threads for each channel.
//! 
//! Work in Progress

extern crate libc;


use std::sync::mpsc;

use crate::acars::common::Reception;

const _SAMPLE_RATE: u32   = 20_000_000;
const _SIZE_RATE_MULTIPLIER: usize = 1; 


// #[repr(C)]
// pub struct hackrf_device;

// #[derive(Debug)]
// pub struct HackRFDevice {
//     ptr: *mut hackrf_device
// }

#[allow(dead_code)]
pub struct Device {
    index: i32,
    central_frequency: f64,
    sample_rate: u32,
    buffer_size: usize,
    // outputs: Vec<mpsc::Sender<Arc<Vec<Complex<f64>>>>>,
}


impl Device {

    
    pub fn new(_index: i32, _frequencies: Vec<f64>, _output: mpsc::Sender<Reception>) {

        let result = match unsafe { hackrf_init() } {
            0     => Ok(()),
            error => Err(format!("Error {}", error)),
        };
        println!("hackrf_init() result: {:?}", result);

        // let mut device: HackRFDevice = unsafe { std::mem::zeroed() };
        // let result = match unsafe { hackrf_open(&mut device.ptr) } {
        //     0     => Ok(device),
        //     error => Err(format!("Error {}", error)),
        // };
        // println!("hackrf_open() result: {:?}", result);

    }

}


#[link(name = "hackrf")]
extern "C" {
    pub fn hackrf_init() -> libc::c_int;


    // pub fn hackrf_open(device: *mut *mut hackrf_device) -> libc::c_int;

}
