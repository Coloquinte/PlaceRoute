![Build](https://github.com/Coloquinte/PlaceRoute/actions/workflows/build.yml/badge.svg)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/8cfe5dc06da74f399fc007e69b742cdc)](https://www.codacy.com/gh/Coloquinte/PlaceRoute/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=Coloquinte/PlaceRoute&amp;utm_campaign=Badge_Grade)
[![GitHub](https://img.shields.io/github/license/coloquinte/torchsr?color=blue)](https://opensource.org/licenses/MIT)

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

For the latest version of other OSes, install the dependencies and build the module. For example:
``` bash
sudo apt-get install cmake libboost-all-dev libeigen3-dev
cd pycoloquinte
python setup.py install
```

### C++ library

Install dependencies and build Coloquinte using CMake :
``` bash
sudo apt-get install libboost-all-dev libeigen3-dev
cmake -B build; cmake --build build
```

