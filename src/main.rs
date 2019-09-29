
use std::thread;

mod rtl;
mod acars;
mod output;

fn main() {

    let frequencies = vec!(
        129.125,
        // 131.550,
        // 131.725,
    );

    let handle = thread::spawn(|| { rtl::Device::new(0, frequencies); });
    handle.join().unwrap();


}
