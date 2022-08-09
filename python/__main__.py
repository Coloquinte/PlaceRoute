
import argparse

import numpy as np

from .circuit import Circuit

parser = argparse.ArgumentParser()
parser.add_argument("instance", help="Benchmark instance")
parser.add_argument("solution_file", help="Placement result")
args = parser.parse_args()

circuit = Circuit.read_ispd(args.instance)
circuit.place()
circuit.write_pl(args.solution_file)
#circuit.benchmark_quadratic_models()

