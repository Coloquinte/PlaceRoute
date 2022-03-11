
import sys
import argparse

import numpy as np

from .ispd_reader import read_ispd


parser = argparse.ArgumentParser(description='Benchmark quadratic net models.')
parser.add_argument('circuit', type=str)
parser.add_argument('--model', choices=['STAR', 'BSTAR', 'B2B'], required=True)
parser.add_argument('--nb-steps', type=int, default=20)
parser.add_argument('--epsilon', type=float, default=1.0)
parser.add_argument('--relaxation', type=float, default=0.0)

args = parser.parse_args()
circuit = read_ispd(args.circuit)
circuit.benchmark_quadratic_models(args.model, args.nb_steps, args.epsilon, args.relaxation)

