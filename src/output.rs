
use std::sync::mpsc;
use crate::chrono::Timelike;
use crate::acars::common::Reception;

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

mod cli {

}

mod network {

}