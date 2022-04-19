
import sys

import numpy as np

from .ispd_reader import read_ispd

circuit = read_ispd(sys.argv[1])
circuit.place()
#circuit.benchmark_quadratic_models()

