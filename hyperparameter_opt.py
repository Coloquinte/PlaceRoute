import argparse
import math
import os
import time

import coloquinte

import numpy as np

from openbox import space as sp
from openbox import Optimizer


class BenchmarkRun:
    def __init__(
        self,
        benchmark,
        global_params=None,
        detailed_params=None,
        ignore_macros=False,
        prefix="benchmarks/ISPD06/",
    ):
        if global_params is None:
            global_params = coloquinte.GlobalPlacerParameters()
        if detailed_params is None:
            detailed_params = coloquinte.DetailedPlacerParameters()
        self._global_params = global_params
        self._detailed_params = detailed_params
        self._benchmark = benchmark
        self._ignore_macros = ignore_macros
        self._prefix = prefix

    @property
    def gp(self):
        return self._global_params

    @property
    def dp(self):
        return self._detailed_params

    @property
    def benchmark(self):
        return self._benchmark

    @property
    def ignore_macros(self):
        return self._ignore_macros

    def run(self):
        circuit = coloquinte.Circuit.read_ispd(
            os.path.join(self._prefix, self.benchmark),
            ignore_obstructions=self.ignore_macros,
        )
        # Global placement
        global_start_time = time.time()
        circuit.place_global(self.gp)
        global_end_time = time.time()
        # Detailed placement
        detailed_start_time = time.time()
        circuit.place_detailed(self.dp)
        detailed_end_time = time.time()
        # Return the metrics
        global_duration = global_end_time - global_start_time
        detailed_duration = detailed_end_time - detailed_start_time
        hpwl = circuit.hpwl()
        metrics = {
            "time_total": global_duration + detailed_duration,
            "time_global": global_duration,
            "time_detailed": detailed_duration,
            "hpwl": hpwl,
        }
        return metrics

    @staticmethod
    def _setattr(params, k, v):
        """
        Set a parameter from a string or scalar, including Pybind enums
        """
        old = getattr(params, k)
        if not isinstance(old, (int, float, str)):
            v = type(old).__members__[v]
        setattr(params, k, v)

    @staticmethod
    def from_dict(d):
        gp = coloquinte.GlobalPlacerParameters()
        dp = coloquinte.DetailedPlacerParameters()
        for k, v in d.items():
            if k.startswith("global_"):
                k = k[7:]
                BenchmarkRun._setattr(gp, k, v)
            if k.startswith("detailed_"):
                k = k[9:]
                BenchmarkRun._setattr(dp, k, v)
        return BenchmarkRun(d["benchmark"], gp, dp, d["ignore_macros"], d["prefix"])

    def to_dict(self):
        ret = {}
        for name, params in [("global", self.gp), ("detailed", self.dp)]:
            for p in params.__dir__():
                if p.startswith("__"):
                    continue
                if p == "check":
                    continue
                val = getattr(params, p)
                if not isinstance(val, (int, float, str)):
                    val = val.name
                ret[name + "_" + p] = val
        ret["benchmark"] = self.benchmark
        ret["prefix"] = self._prefix
        ret["ignore_macros"] = self._ignore_macros
        return ret


optimization_variables = [
    sp.Int("global_max_nb_steps", 30, 200, log=True),
    sp.Int("global_nb_initial_steps", 0, 2),
    sp.Real("global_gap_tolerance", 0.02, 0.2, log=True),
    sp.Real("global_distance_tolerance", 2, 10, log=True),
    sp.Real("global_initial_penalty", 0.02, 0.08, log=True),
    sp.Real("global_penalty_update_factor", 1.05, 1.3, log=True),
    sp.Real("global_penalty_cutoff_distance", 2.0, 50.0, log=True),
    sp.Real("global_approximation_distance", 0.1, 10.0, log=True),
    sp.Int("global_max_nb_conjugate_gradient_steps", 100, 1000, log=True),
    sp.Int("global_rough_legalization_nb_steps", 1, 3),
    sp.Int("global_rough_legalization_reopt_square_size", 1, 3),
    sp.Int("global_rough_legalization_reopt_length", 2, 8, log=True),
    sp.Real("global_export_weighting", 0.5, 1.0),
    sp.Int("detailed_nb_passes", 1, 5),
    sp.Int("detailed_local_search_nb_neighbours", 2, 16),
    sp.Int("detailed_local_search_nb_rows", 1, 3),
    sp.Int("detailed_shift_nb_rows", 2, 8),
    sp.Int("detailed_shift_max_nb_cells", 20, 150, log=True),
    sp.Int("detailed_reordering_nb_rows", 1, 3),
    sp.Int("detailed_reordering_max_nb_cells", 1, 8),
]


def show_variables():
    print("Optimization variables:")
    for v in optimization_variables:
        print(f"\t{v.name}")


class HPOptimizer:
    def __init__(
        self,
        args,
        prefix="benchmarks/ISPD06/",
    ):
        benchmarks = args.benchmarks
        if benchmarks is None:
            benchmarks = os.listdir(prefix)
        self.benchmarks = benchmarks

        variables = args.variables
        if variables is None:
            self.variables = optimization_variables
        else:
            self.variables = [v for v in optimization_variables if v.name in variables]

        self.args = args
        self.prefix = prefix

    def default_params():
        data = BenchmarkRun("").to_dict()
        del data["benchmark"]
        del data["ignore_macros"]
        del data["prefix"]
        del data["global_seed"]
        del data["detailed_seed"]
        return data

    def evaluate(self, params_dict):
        """
        Return the geometric mean of the metrics across the benchmarks, with the given time/quality tradeoff
        """
        metrics = {}
        for benchmark in sorted(self.benchmarks):
            for seed in range(self.args.nb_seeds):
                params = dict(params_dict)
                params["benchmark"] = benchmark
                params["ignore_macros"] = self.args.ignore_macros
                params["prefix"] = self.prefix
                params["global.seed"] = seed
                params["detailed.seed"] = seed
                cur = BenchmarkRun.from_dict(params).run()
                for k, v in cur.items():
                    if k in metrics:
                        metrics[k].append(v)
                    else:
                        metrics[k] = [v]
        # Geometric mean
        m = {}
        for k, v in metrics.items():
            m[k] = np.exp(np.mean(np.log(v)))
        return m

    def evaluate_openbox(self, config):
        params_dict = HPOptimizer.default_params()
        for k, v in config.get_dictionary().items():
            params_dict[k] = v
        m = self.evaluate(params_dict)
        t = m["time_total"]
        q = m["hpwl"]
        print("New incumbent evaluated: ")
        for k, v in config.get_dictionary().items():
            print(f"\t{k}: {v}")
        print(f"Objective:\tQuality {q:.0f}\tTime {t:.2f}")
        return {"objs": [q, t]}

    def run(self):
        space = sp.Space()
        space.add_variables(self.variables)
        opt = Optimizer(
            self.evaluate_openbox,
            space,
            num_objs=2,
            num_constraints=0,
            surrogate_type='prf',
            acq_type='ehvi',
            acq_optimizer_type='random_scipy',
            initial_runs=2*len(self.variables)+1,
            max_runs=self.args.max_nb_runs,
            runtime_limit=self.args.time_limit,
            time_limit_per_trial=self.args.time_limit_per_trial,
            task_id='coloquinte_hyperparameter',
            ref_point = [1.0e8, 600],
            random_state=1,
        )
        history = opt.run()
        import pdb; pdb.set_trace()

parser = argparse.ArgumentParser()
parser.add_argument("--show-variables", action="store_true")
parser.add_argument(
    "--benchmarks", help="Benchmark instances to optimize", nargs="+", type=str
)
parser.add_argument(
    "--ignore-macros",
    help="Ignore the macros during placement for more stable benchmarks",
    action="store_true",
)
parser.add_argument(
    "--nb-seeds",
    help="Number of seedsfor each benchmark to reduce noise",
    type=int,
    default=1,
)
parser.add_argument("--max-nb-runs", help="Maximum number of runs", type=int, default=100)
parser.add_argument("--time-limit", help="Time limit for optimization", type=int)
parser.add_argument("--time-limit-per-trial", help="Time limit for each run", type=int, default=600)
parser.add_argument(
    "--variables",
    help="Variables to optimize over",
    type=str,
    nargs="+",
    choices=[v.name for v in optimization_variables],
    metavar="VARIABLES"
)

args = parser.parse_args()
if args.show_variables:
    show_variables()
else:
    optimizer = HPOptimizer(args)
    optimizer.run()
