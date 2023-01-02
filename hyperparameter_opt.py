import argparse
import math
import os
import pickle
import sys
import time

import coloquinte

import numpy as np

import openbox
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
            "metrics.time_total": global_duration + detailed_duration,
            "metrics.time_global": global_duration,
            "metrics.time_detailed": detailed_duration,
            "metrics.hpwl": hpwl,
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
        if "effort" in d:
            effort = d["effort"]
        else:
            effort = 3
        gp = coloquinte.GlobalPlacerParameters(effort)
        dp = coloquinte.DetailedPlacerParameters(effort)
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


def enum_values(enum_type):
    return sorted([v for v in enum_type.__members__.keys()])


optimization_variables = [
    sp.Int("effort", 1, 9),
    sp.Real("global_approximation_distance", 0.5, 10.0, log=True),
    sp.Real("global_approximation_distance_update_factor", 0.95, 1.05, log=True),
    sp.Real("global_conjugate_gradient_error_tolerance", 1.0e-8, 1.0e-4, log=True),
    sp.Real("global_distance_tolerance", 2, 10, log=True),
    sp.Real("global_export_weighting", 0.5, 1.0),
    sp.Real("global_gap_tolerance", 0.02, 0.15, log=True),
    sp.Real("global_initial_penalty", 0.02, 0.08, log=True),
    sp.Int("global_max_nb_conjugate_gradient_steps", 100, 1000, log=True),
    sp.Int("global_max_nb_steps", 20, 500, log=True),
    sp.Int("global_nb_initial_steps", 0, 2),
    sp.Int("global_nb_steps_per_legalization", 1, 2),
    sp.Categorical("global_net_model", enum_values(coloquinte.NetModel)),
    sp.Real("global_penalty_area_exponent", 0.5, 1.0),
    sp.Real("global_penalty_cutoff_distance", 5.0, 100.0, log=True),
    sp.Real("global_penalty_cutoff_distance_update_factor", 0.95, 1.05, log=True),
    sp.Real("global_penalty_update_factor", 1.04, 1.25, log=True),
    sp.Real("global_rough_legalization_bin_size", 3.0, 10.0, log=True),
    sp.Real("global_rough_legalization_coarsening_limit", 0.5, 100.0, log=True),
    sp.Categorical("global_rough_legalization_cost_model", enum_values(coloquinte.LegalizationModel)),
    sp.Int("global_rough_legalization_nb_steps", 1, 5),
    sp.Real("global_rough_legalization_quadratic_penalty", 0.0, 0.005),
    sp.Int("global_rough_legalization_reopt_length", 2, 16, log=True),
    sp.Int("global_rough_legalization_reopt_square_size", 1, 4),
    sp.Real("global_rough_legalization_side_margin", 0.0, 1.5),
    sp.Real("detailed_legalization_ordering_width", 0.0, 1.0),
    sp.Real("detailed_legalization_ordering_y", 0.0, 0.001),
    sp.Int("detailed_local_search_nb_neighbours", 2, 16),
    sp.Int("detailed_local_search_nb_rows", 1, 3),
    sp.Int("detailed_nb_passes", 1, 5),
    sp.Int("detailed_reordering_max_nb_cells", 1, 8),
    sp.Int("detailed_reordering_nb_rows", 1, 3),
    sp.Int("detailed_shift_max_nb_cells", 20, 150, log=True),
    sp.Int("detailed_shift_nb_rows", 2, 8),
]

default_vars = [
    "effort",
    #"global_approximation_distance",
    #"global_conjugate_gradient_error_tolerance",
    #"global_distance_tolerance",
    #"global_export_weighting",
    #"global_gap_tolerance",
    #"global_initial_penalty",
    #"global_max_nb_conjugate_gradient_steps",
    #"global_max_nb_steps",
    #"global_nb_initial_steps",
    #"global_nb_steps_per_legalization",
    #"global_net_model",
    #"global_penalty_area_exponent",
    #"global_penalty_cutoff_distance",
    #"global_penalty_update_factor",
    #"global_rough_legalization_bin_size",
    #"global_rough_legalization_coarsening_limit",
    #"global_rough_legalization_cost_model",
    #"global_rough_legalization_nb_steps",
    #"global_rough_legalization_quadratic_penalty",
    #"global_rough_legalization_reopt_length",
    #"global_rough_legalization_reopt_square_size",
    #"global_rough_legalization_side_margin",
    #"detailed_legalization_ordering_width",
    #"detailed_legalization_ordering_y",
    #"detailed_local_search_nb_neighbours",
    #"detailed_local_search_nb_rows",
    #"detailed_nb_passes",
    #"detailed_reordering_max_nb_cells",
    #"detailed_reordering_nb_rows",
    #"detailed_shift_max_nb_cells",
    #"detailed_shift_nb_rows",
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

    def load_history(self):
        history_file = self.args.history_file
        if os.path.exists(history_file):
            with open(history_file, 'rb') as f:
                return pickle.load(f)
        else:
            return []

    def save_history(self, call_history):
        history_file = self.args.history_file
        with open(history_file, 'wb') as f:
            pickle.dump(call_history, f)

    def is_old_run(self, params, candidate):
        for k, v in params.items():
            if candidate[k] != v:
                return False
        return True

    def find_history(self, params):
        call_history = self.load_history()
        for candidate in call_history:
            if self.is_old_run(params, candidate):
                print(f"Found candidate {candidate}")
                return candidate
        return None

    def history_as_conf(self, config_space):
        history = self.load_history()
        confs = []
        for h in history:
            conf_dict = {}
            for k, v in h.items():
                if k in config_space:
                    conf_dict[k] = v
            conf = openbox.utils.config_space.Configuration(config_space, values=conf_dict)
            confs.append(conf)
        return confs

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
                old = self.find_history(params)
                if old is not None:
                    params = old
                else:
                    cur = BenchmarkRun.from_dict(params).run()
                    for k, v in cur.items():
                        params[k] = v
                    call_history = self.load_history()
                    call_history.append(params)
                    self.save_history(call_history)
                for k in ["metrics.time_total", "metrics.hpwl", "metrics.time_global", "metrics.time_detailed"]:
                    if k in metrics:
                        metrics[k].append(params[k])
                    else:
                        metrics[k] = [params[k]]
        # Geometric mean
        m = {}
        for k, v in metrics.items():
            m[k] = np.exp(np.mean(np.log(v)))
        return m

    def evaluate_openbox(self, config, multiobj):
        print("Evaluating new incumbent:")
        for k, v in config.get_dictionary().items():
            print(f"\t{k}: {v}")
        sys.stdout.flush()
        params_dict = {}
        for k, v in config.get_dictionary().items():
            params_dict[k] = v
        m = self.evaluate(params_dict)
        t = m["metrics.time_total"]
        q = m["metrics.hpwl"]
        sys.stdout.flush()
        if multiobj:
            print(f"Objective:\tQuality {q:.0f}\tTime {t:.2f}")
            return {"objs": [q, t]}
        else:
            factor = np.log(1 + 0.01 * self.args.percents_per_hour) / 3600
            blended = np.exp(np.log(q) + t * factor)
            print(f"Objective:\tQuality {q:.0f}\tTime {t:.2f}\tBlended {blended:.2f}")
            return {"objs": [q, t]}

    def evaluate_multi_objective(self, config):
        return self.evaluate_openbox(config, True)

    def evaluate_single_objective(self, config):
        return self.evaluate_openbox(config, False)

    def run(self):
        space = sp.Space()
        space.add_variables(self.variables)
        initial_configs = self.history_as_conf(space)
        print(f"Found {len(initial_configs)} configurations already evaluated")
        if self.args.percents_per_hour is None:
            opt = Optimizer(
                self.evaluate_multi_objective,
                space,
                num_objs=2,
                num_constraints=0,
                surrogate_type='prf',
                acq_type='ehvi',
                initial_runs=2*len(self.variables)+1,
                max_runs=self.args.max_nb_runs,
                runtime_limit=self.args.time_limit,
                time_limit_per_trial=self.args.time_limit_per_trial,
                task_id='coloquinte_hyperparameter',
                ref_point = [self.args.ref_hpwl, self.args.ref_time],
                random_state=1,
                initial_configurations=initial_configs,
            )
        else:
            opt = Optimizer(
                self.evaluate_single_objective,
                space,
                num_objs=1,
                num_constraints=0,
                surrogate_type='prf',
                initial_runs=2*len(self.variables)+1,
                max_runs=self.args.max_nb_runs,
                runtime_limit=self.args.time_limit,
                time_limit_per_trial=self.args.time_limit_per_trial,
                task_id='coloquinte_hyperparameter',
                random_state=1,
                initial_configurations=initial_configs,
            )
        history = opt.run()
        history.save_json("history.json")

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
parser.add_argument("--percents-per-hour", help="Conversion factor to get one single objective", type=float)
parser.add_argument("--ref-hpwl", help="Reference value for the hpwl in multiobjective", type=float, default=1.0e8)
parser.add_argument("--ref-time", help="Reference value for the time in multiobjective", type=float, default=5000)
parser.add_argument(
    "--variables",
    help="Variables to optimize over",
    type=str,
    nargs="+",
    choices=[v.name for v in optimization_variables],
    metavar="VARIABLES",
    default=default_vars
)
parser.add_argument("--history-file", help="History file to save/load", type=str, default="history.pkl")

args = parser.parse_args()
if args.show_variables:
    show_variables()
else:
    optimizer = HPOptimizer(args)
    optimizer.run()
