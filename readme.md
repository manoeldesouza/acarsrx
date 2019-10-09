# acarsrx
Copyright (c) 2019 Manoel Souza <manoel.desouza@outlook.com.br>

ACARS decoder for RTL-SDR and (in future) other SDR devices

Inspired on Acarsdec Copyright (c) by Thierry Thierry Leconte

https://github.com/TLeconte/acarsdec
and Acars-sdrplay Copyright (c) by Jan van Katwijk <J.vanKatwijk@gmail.com>

https://github.com/JvanKatwijk/acars-sdrplay


 ## Introduction:
acarsrx is designed to be able to run each SDR device in separate thread enabling multiple devices 
to be controlled from a single execution. For each SDR device controlled, multiple channels can be 
used to demodulate and decode multiple frequencies. Special attention is put into decoupling the 
SDR handling and the channels processing (where demodulation, decoding and block formating). This 
way, it is planned for other SDR devices to be controlled from the same application (HackRF is a 
work in progress). Finally, acarsrx is designed in a way to recover from SDR devices disconnects.
Being USB devices the author experienced common disconnect situations and enable the application 
to recover from such situation proved to be an asset. 


## Design:

Each SDR device is planned to have all instructions handling from a single Rust file (as done for 
rtlsdr.rs and being done for hackrf.rs).

Following the same concept, each protocol is planned to be fully defined in a single Rust file (as 
done for acars.rs). acars.rs defines all instructions for plain old ACARS (POA) and is planned to 
also include VDL Mode 2 decoding in future releases.
The system is designed to initiate threads to handle each SDR device and other threads for each 
channel from within the SDR thread. 


### Thread execution example: 

```bash
acarsrx --rtlsdr 0 131.550 131.725 --rtlsdr 1 129.125 129.425
```

```text
                        / acars.rs 131.550 \
          -- rtlsdr.rs 0                    \
        /               \ acars.rs 131.725 \ \
main.rs                                      output.rs
        \               / acars.rs 129.125 / /
          -- rtlsdr.rs 1                    /
                        \ acars.rs 129.425 /
```
### Sample Output:
```text
07:41:55 RADIO-r0-0 131.550 MHz -24 ↗ ⊝ [2........SQ..01XAYULCYUL1ARINC]
07:42:07 RADIO-r0-0 129.125 MHz -21 ↗ ⊝ [2.C-FTJQ1_dY]
07:42:39 RADIO-r0-0 131.725 MHz -14 ↘ ⊝ [2..N8767K_d5.S75AZD8767]
07:45:11 RADIO-r0-0 131.725 MHz -23 ↗ ⊝ [2........SQ..02XSYULCYUL04527N07344WV136975/]
07:46:05 RADIO-r0-0 131.550 MHz -21 ↗ ⊝ [2........SQ..02XAYULCYUL14528N07344WV136975/ARINC]
07:46:36 RADIO-r0-0 131.725 MHz -20 ↘ ⊝ [2..N8767N_d1.S78AZD8767]
-------- ---------- ----------- --- - - -------------------------
    |       |          |         |  | | ACARS Block
    |       |          |         |  | SingleBlock or MultiBlock
    |       |          |         |  Uplink or Downlink
    |       |          |         Signal Strenght
    |       |          Channel Frequency
    |       Station, device & Channel number 
    Time (hh:mm:ss)
``` 


## Challenges:

There is some odd behavior still under investigation by the author. The simple sample_size achived by 
diving the device SAMPLE_RATE per the CHANNEL_RATE, shows poor decoding performance, with improvements 
being a achived with multiple of these values. To achieve better performance the author applied a 
factor constant named SIZE_RATE_MULTIPLIER = 20. Values from 10 in general are enough for good results.
When compared with other acars decoders this seems to be a deffect still to be investigated.


## Next steps:

There are some functionalities still not present in acarsrx to achieve V1.0, like:
1. Multiple Reception qualities to be submitted to output.rs (Bad CRC block, Bad parity bit, etc). It 
   will be for the output.rs thread to define what is to be printed and how.

2. CRC error checking, correction and re parity check.

Other improvements, new functionalities are planned for future:

3. Implement multiple output formats in addition to current single-line print, for example JSON, Ncurses
   based interface, UDP transmission, Save to file. 

4. Enable ACARS Block assembly into ACARS Messages (Type-B messages) with translation of Labels into 
   SMI, and standard message formating (Following AEEC 620-8 specification).

5. Control of HackRF SDR. Other SDR devices can be controled as well given that similar methods to
   enable the IQ buffer to be retrieved and transformed into vectors of complex numbers.

6. Implementation of VDL Mode 2 decoding. By design a single SDR should enable decoding of multiple 
   protocols and both Plain Old ACARS and VDL Mode 2 should be possible to be decoded, give that enough
   bandwidth and channel separation is provided.


## Future developments:

The exercise on implementing a decoupling between ACARS demodulation and decoding from the sample
extraction from SDR device, inspires the author to consider developing a general SDR driver which
would be able to control and generate samples for multiple protocols (couting on community support)
and outputs. Such server would enable multiple decoders to consume the same sources and multiple 
outputs to be triggered from those.

| Input |  Process  |  Output |
|:-----:|:---------:|:-------:|
|RTLSDR | AM FM SSB | speaker |
|HACKRF | ACARS     | wav     |
|AIRSPY | CW        | json    |
|etc    | APRS      | txt     |
|       | etc       | tcp/udp |

Suggested name: sdrD (software defined radio Daemon).

## Additional details on acarsrx:

The modules section below describes the internal functions in acarsrx

My regards, Manoel Souza
08-Oct-2019


ToDo:

 - [ ] Improve Reception object handling in output;
 - [ ] ACARS Frame error checking & correction;


