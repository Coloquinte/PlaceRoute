![Build](https://github.com/Coloquinte/PlaceRoute/actions/workflows/build.yml/badge.svg)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/8cfe5dc06da74f399fc007e69b742cdc)](https://www.codacy.com/gh/Coloquinte/PlaceRoute/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=Coloquinte/PlaceRoute&amp;utm_campaign=Badge_Grade)
[![GitHub](https://img.shields.io/github/license/coloquinte/PlaceRoute?color=blue)](https://opensource.org/licenses/MIT)
[![PyPI](https://img.shields.io/pypi/v/coloquinte?color=orange)](https://pypi.org/project/coloquinte/)

# Coloquinte Place&Route

Coloquinte is a Place&Route tool for electronic circuits.
Its goal is to provide a single package of well-tested and well-tuned Place&Route algorithms, to be used in open source electronic design toolchains.

This replaces and extends the [placement library](https://github.com/Coloquinte/Coloquinte_placement) used in the [Coriolis](https://gitlab.lip6.fr/vlsi-eda/coriolis/) toolchain.

## Using Coloquinte

ISPD placement and routing benchmarks are available directly in this repository using Git LFS. [Placement benchmarks](https://github.com/Coloquinte/PlaceRouteBenchmarks) can be run with the Python package:
``` bash
coloquinte ISPD06/adaptec1
```

For other applications, you can use Coloquinte as a C++ library or a Python package: see `src/coloquinte.hpp` or `help(coloquinte.Circuit)`.

## Installing Coloquinte

### Python package

On Linux, you may install Coloquinte from pip:
``` bash
pip install coloquinte
```

For the latest version or for other OSes, install the dependencies and build the module. For example:
``` bash
sudo apt-get install g++ cmake libboost-all-dev libeigen3-dev liblemon-dev
cd pycoloquinte
python setup.py install
```

### C++ library

Install dependencies and build Coloquinte using CMake:
``` bash
sudo apt-get install g++ cmake libboost-all-dev libeigen3-dev liblemon-dev
cmake -B build; cmake --build build; ctest --test-dir build
```

Or using Meson:
```
meson setup build; meson compile -C build; meson test -C build
```

## Benchmarks

Coloquinte is tested on the [ISPD06 benchmark suite](https://dl.acm.org/doi/10.1145/1123008.1123042). Below is the reported half-perimeter wirelength on these benchmarks (x10<sup>7</sup>) for various effort parameters.
Higher effort = higher quality but higher runtime.


| Benchmark | Effort 1 WL | Effort 3 WL | Effort 6 WL |
| --------- | ----------- | ----------- | ----------- |
| adaptec1  |        8,04 |        7,62 |        7,55 |
| adaptec2  |        8,83 |        8,48 |        8,38 |
| adaptec3  |       21,07 |       20,31 |       19,95 |
| adaptec4  |       18,84 |       17,97 |       17,77 |
| adaptec5  |       32,65 |       31,18 |       30,89 |
| bigblue1  |       10,01 |        9,36 |        9,20 |
| bigblue2  |       14,79 |       14,25 |       14,17 |
| bigblue3  |       35,23 |       32,67 |       31,98 |
| bigblue4  |       83,15 |       79,16 |       77,66 |
| newblue1  |       60,68 |       18,75 |       34,54 |
| newblue2  |       19,38 |       18,31 |       17,97 |
| newblue3  |       33,70 |       26,06 |       25,84 |
| newblue4  |       24,32 |       23,35 |       23,33 |
| newblue5  |       42,29 |       39,39 |       39,02 |
| newblue6  |       48,43 |       45,80 |       45,39 |
| newblue7  |      102,11 |       97,81 |       96,47 |
