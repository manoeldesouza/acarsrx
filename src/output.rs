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

//! Provides all output functions for ACARS blocks. Currently only prints a single-line block
//! representation to the screen, but will in future enable multiple output formats like JSON,
//! detailed Reception and Block information and also forward the Reception detail via UDP to
//! other hosts. Potential future development also include ncurses based interface and / or 
//! HTTP portal.


use std::sync::mpsc;
use crate::chrono::Timelike;
use crate::acars::common::Reception;

/// Initiates the thread controlling the all outputs of the program (only stdout for the time being)
pub fn thread(input: mpsc::Receiver<Reception>) {

    loop {

        let reception = input.recv().expect("Error on OUTPUT");

        let printable_block = match reception.acars_block {
            Some(block) => block.get_essential(),
            None => String::from(""),
        };

        println!("{:0>2}:{:0>2}:{:0>2} {}[{}] {:.3} MHz {: >3} {}", 
            reception.date_time.hour(),
            reception.date_time.minute(),
            reception.date_time.second(),
            reception.station, reception.channel, reception.frequency, reception.level, 
            printable_block);
    }
}
