
![build](https://github.com/Coloquinte/PlaceRoute/actions/workflows/build.yml/badge.svg)
[![GitHub](https://img.shields.io/github/license/coloquinte/torchsr?color=blue)](https://opensource.org/licenses/MIT)

# Coloquinte Place&Route

Coloquinte is a Place&Route tool for electronic circuits.
Its goal is to provide a single package of well-tested and well-tuned Place&Route algorithms, to be used in open source electronic design toolchains.

## Using Coloquinte

Build Coloquinte using CMake :
```
cmake -B build; cmake --build build
```


ISPD placement and routing benchmarks are available directly in this repository using Git LFS. Currently, only placement benchmarks can be run.
```
python -m pycoloquinte benchmarks/ISPD06/adaptec1/adaptec1.aux
```


To use Coloquinte as a library, setup a `coloquinte::Circuit` object, and call `coloquinte::place`.
