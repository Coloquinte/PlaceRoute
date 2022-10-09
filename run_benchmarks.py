import argparse
import math
import os
import sqlite3
import time

import coloquinte

import numpy as np


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


class BlackboxFloatVariable:
    def __init__(self, min, max):
        self.min = min
        self.max = max

    def define(self, model):
        return model.float(self.min, self.max)

    def value(self, model_value):
        return model_value

    def model_value(self, value):
        return value


class BlackboxLogFloatVariable:
    def __init__(self, min, max):
        self.min = min
        self.max = max

    def define(self, model):
        return model.float(math.log(self.min), math.log(self.max))

    def value(self, model_value):
        return math.exp(model_value)

    def model_value(self, value):
        return math.log(value)


class BlackboxIntVariable:
    def __init__(self, min, max):
        self.min = min
        self.max = max

    def define(self, model):
        return model.int(self.min, self.max)

    def value(self, model_value):
        return model_value

    def model_value(self, value):
        return value


class BlackboxLogIntVariable:
    def __init__(self, min, max):
        self.min = min
        self.max = max

    def define(self, model):
        return model.float(math.log(self.min), math.log(self.max))

    def value(self, model_value):
        return int(round(math.exp(model_value)))

    def model_value(self, value):
        return math.log(value)


class BlackboxEnumVariable:
    def __init__(self, enum):
        self.enum = enum

    def define(self, model):
        return model.int(0, len(self.enum.__members__))

    def value(self, model_value):
        # TODO
        return self.enum(model_value)

    def model_value(self, value):
        # TODO
        return value


optimization_variables = {
    "global_max_nb_steps": BlackboxIntVariable(20, 150),
    "global_gap_tolerance": BlackboxLogFloatVariable(0.01, 0.2),
    "global_initial_penalty": BlackboxFloatVariable(0.01, 0.05),
    "global_penalty_update_factor": BlackboxLogFloatVariable(1.01, 1.3),
    "global_penalty_cutoff_distance": BlackboxLogFloatVariable(2.0, 50.0),
    "global_approximation_distance": BlackboxLogFloatVariable(0.1, 10.0),
    "global_max_nb_conjugate_gradient_steps": BlackboxLogIntVariable(100, 1000),
    "global_nb_rough_legalization_steps": BlackboxIntVariable(1, 3),
    "global_export_weighting": BlackboxFloatVariable(0.5, 1.0),
    "detailed_nb_passes": BlackboxIntVariable(1, 3),
    "detailed_local_search_nb_neighbours": BlackboxIntVariable(1, 6),
    "detailed_local_search_nb_rows": BlackboxIntVariable(1, 3),
    "detailed_shift_nb_rows": BlackboxIntVariable(2, 10),
    "detailed_shift_max_nb_cells": BlackboxLogIntVariable(20, 200),
}


class Optimizer:
    def __init__(
        self,
        variables,
        time_for_quality,
        benchmarks=None,
        nb_runs=1,
        ignore_macros=False,
        prefix="benchmarks/ISPD06/",
    ):
        if benchmarks is None:
            benchmarks = os.listdir(prefix)
        if variables is None:
            variables = list(optimization_variables.keys())
        self.variable_names = variables
        self.benchmarks = benchmarks
        self.nb_runs = nb_runs
        self.ignore_macros = ignore_macros
        self.prefix = prefix
        self.time_for_quality = time_for_quality

    def default_params():
        data = BenchmarkRun("").to_dict()
        del data["benchmark"]
        del data["ignore_macros"]
        del data["prefix"]
        del data["global_seed"]
        del data["detailed_seed"]
        return data

    def evaluate_single(self, params_dict):
        return BenchmarkRun.from_dict(params_dict).run()

    def evaluate(self, params_dict):
        """
        Return the geometric mean of the metrics across the benchmarks, with the given time/quality tradeoff
        """
        metrics = {}
        for benchmark in sorted(self.benchmarks):
            for seed in range(self.nb_runs):
                params = dict(params_dict)
                params["benchmark"] = benchmark
                params["ignore_macros"] = self.ignore_macros
                params["prefix"] = self.prefix
                params["global.seed"] = seed
                params["detailed.seed"] = seed
                cur = self.evaluate_single(params)
                for k, v in cur.items():
                    if k in metrics:
                        metrics[k].append(v)
                    else:
                        metrics[k] = [v]
        # Geometric mean
        m = {}
        for k, v in metrics.items():
            m[k] = np.exp(np.mean(np.log(v)))
        t = m["time_total"]
        q = m["hpwl"]
        # Just normalization so numbers are not too horrible
        norm_q = 1.0e6
        norm_t = 100
        return math.exp(
            math.log(q / norm_q) + math.log(t / norm_t) / self.time_for_quality
        )

    def evaluate_localsolver(self, params_map):
        params_dict = Optimizer.default_params()
        print("Evaluation of a new incumbent: ")
        for i, name in enumerate(self.variable_names):
            v = optimization_variables[name]
            params_dict[name] = v.value(params_map[i])
            print(f"\t{name}: {params_dict[name]}")
        ret = self.evaluate(params_dict)
        print(f"Objective function: {ret}")
        return ret

    def define_variables(self, model):
        return [
            optimization_variables[name].define(model) for name in self.variable_names
        ]

    def run(self, time_limit=None):
        import localsolver

        with localsolver.LocalSolver() as ls:
            # Create a simple model with an external function
            model = ls.model
            model_params = self.define_variables(model)
            f = model.create_double_external_function(self.evaluate_localsolver)
            func_call = model.call(f, *model_params)
            model.minimize(func_call)
            surrogate_params = f.external_context.enable_surrogate_modeling()
            model.close()

            # Setup the initial values as the default parameters
            params_dict = Optimizer.default_params()
            for name, p in zip(self.variable_names, model_params):
                v = optimization_variables[name]
                p.value = v.model_value(params_dict[name])

            if time_limit is not None:
                ls.param.time_limit = time_limit

            # Solve
            ls.solve()

            # Get the result
            result = {}
            for name, p in zip(self.variable_names, model_params):
                v = optimization_variables[name]
                result[name] = v.value(p.value)
            print(f"Optimization result: {result}")


parser = argparse.ArgumentParser()
parser.add_argument(
    "benchmarks", help="Benchmark instances to optimize", nargs="+", type=str
)
parser.add_argument(
    "--ignore-macros",
    help="Ignore the macros during placement for more stable benchmarks",
    action="store_true",
)
parser.add_argument(
    "--nb-runs",
    help="Number of runs for each benchmark to reduce noise",
    type=int,
    default=1,
)
parser.add_argument(
    "--time-quality-tradeoff",
    help="Tradeoff between time and quality; higher is higher quality",
    type=float,
    default=10.0,
)
parser.add_argument("--time-limit", help="Time limit for optimization", type=int)
parser.add_argument(
    "--variables",
    help="Variables to optimize over",
    type=str,
    nargs="+",
    choices=[v for v in optimization_variables.keys()],
)

args = parser.parse_args()

optimizer = Optimizer(
    variables=args.variables,
    time_for_quality=args.time_quality_tradeoff,
    benchmarks=args.benchmarks,
    nb_runs=args.nb_runs,
    ignore_macros=args.ignore_macros,
)
optimizer.run(time_limit=args.time_limit)
