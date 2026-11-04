# Plasma Bubble Indices for the NovAtel OEM719

Three small C programs that measure ionospheric disturbances from raw GNSS receiver data.

Equatorial plasma bubbles are regions of depleted plasma in the night-time ionosphere. They disturb the GNSS signals passing through them, which makes them visible to an ordinary receiver. These programs turn logged OEM719 measurements into the three standard indices used to detect that disturbance.

## The programs

| Program | Index | Measurement used | What it shows |
|---|---|---|---|
| `S4.c` | S4 | C/N0 | Amplitude scintillation |
| `adr.c` | σφ | L1 carrier phase | Phase scintillation |
| `roti.c` | ROTI | L1 and L2 carrier phase | TEC irregularities, the bubble itself |

Scintillation is a side effect of a bubble and is strongest near its edges. ROTI follows the electron content directly, but it needs a receiver tracking two frequencies.

## Quick start

Each program is a single file that needs only `gcc`.

```bash
gcc -O2 -o s4 S4.c -lm
gcc -O2 -o adr adr.c -lm
gcc -O2 -o roti roti.c -lm
```

Pass your CSV file to any of them, and results print to the terminal.

```bash
./s4 data.csv
./adr data.csv
./roti data.csv
```

With no argument, each program looks for a file named `lab_5hz.csv`.

## Input format

The programs read a CSV with one row per satellite signal per epoch, in time order. The first line is a header and is skipped. The columns are read by position, and any columns after these 12 are ignored.

| # | Column | Meaning | S4 | adr | roti |
|---|---|---|:-:|:-:|:-:|
| 1 | `time_s` | GPS seconds of week | ✓ | ✓ | ✓ |
| 2 | `week` | GPS week, whole number | – | – | – |
| 3 | `system` | 0 = GPS, 3 = Galileo, 5 = QZSS | ✓ | ✓ | ✓ |
| 4 | `prn` | Satellite number | ✓ | ✓ | ✓ |
| 5 | `signal` | Signal type, see below | ✓ | ✓ | ✓ |
| 6 | `cn0_dbhz` | C/N0 in dB-Hz | ✓ | ✓ | – |
| 7 | `adr_cycles` | Carrier phase in cycles | – | ✓ | ✓ |
| 8 | `doppler_hz` | Doppler in Hz | – | – | – |
| 9 | `psr_m` | Pseudorange in metres | – | – | – |
| 10 | `locktime_s` | Lock time in seconds | – | ✓ | ✓ |
| 11 | `phase_lock` | 1 = phase locked | ✓ | ✓ | ✓ |
| 12 | `parity_known` | 1 = parity known | – | ✓ | ✓ |

A dash means the program doesn't use the value, but the column must still be there.

The signal codes come from the channel tracking status word in NovAtel's `RANGE` log:

| System | L1 signal | L2 signal |
|---|---|---|
| GPS | 0 (L1 C/A) | 17 (L2C) |
| Galileo | 2 (E1) | – |
| QZSS | 0 (L1 C/A) | – |

The parser that makes this CSV from binary `RANGEB` and `RANGECMPB` logs is not part of this repo, and neither is the lab data.

## How each index is calculated

### S4

S4 is the normalized standard deviation of signal power, calculated for every satellite and signal over 60 s windows (Kintner et al., 2007).

```text
S4 = √( (⟨I²⟩ − ⟨I⟩²) / ⟨I⟩² )
I  = 10^( (C/N0) / 10 )
```

The levels follow González (2020):

| S4 | Printed label |
|---|---|
| 0.1 or below | No label |
| Above 0.1 | `weak` |
| Above 0.5 | `STRONG` |

### σφ

σφ is the standard deviation of the detrended carrier phase in radians, calculated for every L1 satellite over 60 s windows (Kintner et al., 2007).

1. The satellite with the highest mean C/N0 becomes the reference, and subtracting it removes the receiver clock.
2. A straight line fitted over 10 s removes the slow trend, which acts as a high-pass filter at roughly 0.1 Hz. A line is used instead of a moving average because the samples are unevenly spaced.
3. Windows with a gap over 1 s or a cycle slip are printed as skipped.

σφ has no labelled levels because thresholds for this setup haven't been chosen yet.

### ROTI

ROTI is the standard deviation of the rate of change of TEC, calculated for every GPS satellite with both L1 and L2 (Pi et al., 1997).

1. TEC comes from the L1 and L2 carrier phases, using samples on whole 30 s epochs. This is the sampling used in the papers, and ROTI changes strongly with sample rate (Jacobsen, 2014).
2. ROT is the change in TEC per minute between samples exactly 30 s apart with no cycle slip.
3. ROTI is calculated over 5 min windows.

```text
ROTI = √( ⟨ROT²⟩ − ⟨ROT⟩² )
```

The levels follow González (2020):

| ROTI (TECU/min) | Printed label |
|---|---|
| Below 0.25 | `none` |
| 0.25 to below 0.5 | `weak` |
| 0.5 to below 1 | `moderate` |
| 1 or above | `STRONG` |

All three programs drop the last unfinished window.

## Example output

These S4 and σφ results come from a real lab log recorded in Calgary.

```text
$ ./s4 data.csv
27900 rows from data.csv
sys prn sig   start_s    samples   S4
  0  31   0  591346.8     300    0.088
  0  31   0  591406.8     300    0.133 weak
```

```text
$ ./adr data.csv
15500 L1 rows from data.csv
reference: sys 3 prn 3, mean C/N0 52.3 dB-Hz
sys prn   start_s    samples   sigma_phi (rad)
  0  31  591346.8     300    0.070
  0  31  591406.8     300    0.068
```

These ROTI lines come from made-up L1 + L2 test data, because the lab receiver only tracked L1.

```text
prn   start_s    ROTs   ROTI (TECU/min)
  3  600030.0    10    0.268 weak
  4  600030.0    10    0.682 moderate
  4  605250.0    10    1.066 STRONG
```

## References

- González, G. de L. (2020). Spread-F occurrence during geomagnetic storms near the southern crest of the EIA in Argentina. Accepted for publication in *Advances in Space Research*.
- Jacobsen, K. S. (2014). The impact of different sampling rates and calculation time intervals on ROTI values. *Journal of Space Weather and Space Climate*, 4, A33. [doi:10.1051/swsc/2014031](https://doi.org/10.1051/swsc/2014031)
- Kintner, P. M., Ledvina, B. M., and de Paula, E. R. (2007). GPS and ionospheric scintillations. *Space Weather*, 5, S09003. [doi:10.1029/2006SW000260](https://doi.org/10.1029/2006SW000260)
- Pi, X., Mannucci, A. J., Lindqwister, U. J., and Ho, C. M. (1997). Monitoring of global ionospheric irregularities using the worldwide GPS network. *Geophysical Research Letters*, 24(18), 2283–2286.
