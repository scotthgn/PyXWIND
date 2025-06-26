# PyXWIND

Python wrapper for XWIND. Calculates emission line profiles for accretion disc winds.
Uses a semi-analytic approach, hence this is much faster than simulating winds

## Requirements

- Python 3.10 or above (tested Python 3.12.4)
- numpy (tested on version 1.26.4)
- matplotlib (tested on version 3.9.1)
- meson (if using Python 12 or above) - found at: https://mesonbuild.com/
  - If Python 11 or lower, any fortran complior should do
  - (The meson requirement comes from the numpy f2py bindings used to import the XWIND code into python)


## Installation

Currently:
1. Clone the repository
2. From within the main directory run the compile script compile.sh

```
./compile.sh
```

## Documentation

Coming later...
