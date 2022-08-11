
import argparse

import numpy as np

from .circuit import Circuit

parser = argparse.ArgumentParser()
parser.add_argument("instance", help="Benchmark instance")
parser.add_argument("--solution", help="Placement result")
parser.add_argument("--effort", help="Placement effort", type=int, default=3)
args = parser.parse_args()

circuit = Circuit.read_ispd(args.instance)
circuit.place(args.effort)
circuit.write_pl(args.solution)
