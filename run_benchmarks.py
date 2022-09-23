import argparse
import math
import os
import sqlite3
import time

import coloquinte

import numpy as np


class BenchmarkRun:
    def __init__(self, benchmark, global_params=None, detailed_params=None, prefix="benchmarks/ISPD06/"):
        if global_params is None:
            global_params = coloquinte.GlobalPlacerParameters()
        if detailed_params is None:
            detailed_params = coloquinte.DetailedPlacerParameters()
        self._global_params = global_params
        self._detailed_params = detailed_params
        self._benchmark = benchmark
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

    def run(self):
        circuit = coloquinte.Circuit.read_ispd(
            os.path.join(self._prefix, self.benchmark))
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
            "hpwl": hpwl
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
        return BenchmarkRun(d["benchmark"], gp, dp, d["prefix"])

    def to_dict(self):
        ret = {}
        for name, params in [("global", self.gp), ("detailed", self.dp)]:
            for p in params.__dir__():
                if p.startswith("__"):
                    continue
                val = getattr(params, p)
                if not isinstance(val, (int, float, str)):
                    val = val.name
                ret[name + "_" + p] = val
        ret["benchmark"] = self.benchmark
        ret["prefix"] = self._prefix
        return ret


class BenchmarkRunner:
    def __init__(self, benchmarks=None, time_for_quality=None, prefix="benchmarks/ISPD06/"):
        if benchmarks is None:
            benchmarks = os.listdir(prefix)
        self.benchmarks = benchmarks
        self.prefix = prefix
        self.time_for_quality = time_for_quality

    @staticmethod
    def default_params():
        data = BenchmarkRun("").to_dict()
        del data["benchmark"]
        del data["prefix"]
        return data

    @staticmethod
    def columns():
        data = BenchmarkRun("").to_dict()
        columns = list(data.keys())
        return BenchmarkRunner.metric_columns() + columns

    @staticmethod
    def metric_columns():
        return ["time_total", "time_global", "time_detailed", "hpwl"]

    @staticmethod
    def param_columns():
        data = BenchmarkRun("").to_dict()
        param_columns = list(k for k in data.keys() if k.startswith(
            "global") or k.startswith("detailed"))
        return param_columns

    @staticmethod
    def unique_columns():
        return ["benchmark", "prefix"] + BenchmarkRunner.param_columns()

    @staticmethod
    def create_db():
        # It's ugly from an SQL standpoint, but this way the table includes any field we might need
        columns = ', '.join(BenchmarkRunner.columns())
        unique_columns = ', '.join(BenchmarkRunner.unique_columns())
        sql = f"CREATE TABLE IF NOT EXISTS ColoquinteBenchmarks ({columns}, UNIQUE ({unique_columns}))"
        con = sqlite3.connect("benchmarks.db")
        c = con.cursor()
        c.execute(sql)
        con.commit()
        c.close()

    def get_all_params(self):
        """
        Return all placement parameters that have been used
        """
        param_columns = BenchmarkRunner.param_columns()
        param_columns_sql = ", ".join(param_columns)
        sql = f"SELECT DISTINCT {param_columns_sql} FROM ColoquinteBenchmarks"
        con = sqlite3.connect("benchmarks.db")
        c = con.cursor()
        c.execute(sql)
        ret = c.fetchall()
        return [{k: v for k, v in zip(param_columns, r)} for r in ret]

    def get_data(self, data):
        """
        Get the optimization metrics from the parameters
        """
        key_columns = BenchmarkRunner.unique_columns()
        if sorted(data.keys()) != sorted(key_columns):
            raise RuntimeError("Queried data does not match unique columns")
        values = [data[k] for k in key_columns]
        key_columns_sql = " and ".join([c + "=?" for c in key_columns])
        metric_columns = BenchmarkRunner.metric_columns()
        metric_columns_sql = ", ".join(metric_columns)
        sql = f"SELECT {metric_columns_sql} FROM ColoquinteBenchmarks WHERE {key_columns_sql}"
        con = sqlite3.connect("benchmarks.db")
        c = con.cursor()
        c.execute(sql, values)
        ret = c.fetchall()
        if len(ret) == 0:
            return None
        elif len(ret) == 1:
            assert len(ret[0]) == len(metric_columns)
            return {k: v for k, v in zip(metric_columns, ret[0])}
        else:
            raise RuntimeError("Multiple lines found")

    def save_data(self, data):
        """
        Save the optimization metrics for these parameters
        """
        expected_columns = sorted(BenchmarkRunner.columns())
        all_columns = sorted(data.keys())
        if all_columns != expected_columns:
            raise RuntimeError("Saved data doesn't have expected columns")
        columns = []
        values = []
        placeholders = []
        for k, v in data.items():
            columns.append(k)
            values.append(v)
            placeholders.append("?")
        columns = ', '.join(columns)
        placeholders = ', '.join(placeholders)
        sql = f"INSERT INTO ColoquinteBenchmarks ({columns}) VALUES ({placeholders})"
        con = sqlite3.connect("benchmarks.db")
        c = con.cursor()
        c.execute(sql, values)
        con.commit()
        c.close()

    def run(self, params_dict):
        data = self.get_data(params_dict)
        if data is not None:
            return data
        bench_run = BenchmarkRun.from_dict(params_dict)
        metrics = bench_run.run()
        data = dict(params_dict)
        data.update(metrics)
        self.save_data(data)
        return self.get_data(params_dict)

    def run_metrics(self, params_dict):
        """
        Return a metrics of the run, with a tradeoff between time and quality
        """
        assert self.time_for_quality > 0
        m = self.run(params_dict)
        t = m["time_total"]
        q = m["hpwl"]
        print(
            f"Evaluated benchmark {params_dict['benchmark']} at {q} after {t:.0f}s")
        # Just normalization so numbers are not too horrible
        norm_q = 1.0e6
        norm_t = 100
        return math.exp(math.log(q / norm_q) + math.log(t / norm_t) / self.time_for_quality)

    def run_metrics_all(self, params_dict):
        """
        Return the geometric mean of the metrics across the benchmarks, with the given time/quality tradeoff
        """
        metrics = []
        for benchmark in sorted(self.benchmarks):
            params = dict(params_dict)
            params["benchmark"] = benchmark
            params["prefix"] = self.prefix
            metrics.append(self.run_metrics(params))
        # Geometric mean
        return np.exp(np.mean(np.log(metrics)))

    def evaluate_localsolver(self, params_map):
        params_dict = BenchmarkRunner.default_params()
        param_names = self.get_localsolver_param_names()
        for i in range(len(params_map)):
            params_dict[param_names[i]] = params_map[i]
        print("Evaluation of a new incumbent: ")
        for name in param_names:
            print(f"\t{name}: {params_dict[name]}")
        ret = self.run_metrics_all(params_dict)
        print(f"Objective function: {ret}")
        return ret

    def get_localsolver_param_names(self):
        return [
            "global_max_nb_steps",
            "global_gap_tolerance",
            "global_initial_penalty",
            "global_penalty_update_factor",
            "global_penalty_cutoff_distance",
            "global_approximation_distance",
            "global_max_nb_conjugate_gradient_steps",
            "global_nb_rough_legalization_steps",
            "detailed_nb_passes",
            "detailed_local_search_nb_neighbours",
            "detailed_local_search_nb_rows",
            "detailed_shift_nb_rows",
        ]

    def make_localsolver_parameters(self, model):
        return [
            model.int(20, 60),
            model.float(0.01, 0.2),
            model.float(0.01, 0.05),
            model.float(1.01, 1.3),
            model.float(2.0, 50.0),
            model.float(0.1, 10.0),
            model.int(100, 1000),
            model.int(1, 3),
            model.int(1, 3),
            model.int(1, 6),
            model.int(1, 3),
            model.int(2, 10),
        ]

    def save_localsolver_solutions(self, surrogate_params):
        """
        Register all pre-existing solutions for LocalSolver
        """
        existing = self.get_all_params()
        print(f"{len(existing)} different configurations found")
        nb_added = 0
        for params_dict in existing:
            if self.save_localsolver_solution(surrogate_params, params_dict):
                nb_added += 1
        print(f"{nb_added} configurations added")

    def save_localsolver_solution(self, surrogate_params, params_dict):
        """
        Register an existing solution for LocalSolver
        """
        value = self.run_metrics_all(params_dict)
        params_names = self.get_localsolver_param_names()
        # Do not register if any parameter is not optimized but has a non-default value
        default_dict = BenchmarkRunner.default_params()
        for k, v in params_dict.items():
            if k not in params_names and params_dict[k] != default_dict[k]:
                return False
        # Do not register if any of the benchmarks has not been run
        for benchmark in self.benchmarks:
            benchmark_dict = dict(params_dict)
            benchmark_dict["benchmark"] = benchmark
            benchmark_dict["prefix"] = self.prefix
            if self.get_data(benchmark_dict) is None:
                return False
        # Register the solution
        evaluation_point = surrogate_params.create_evaluation_point()
        for name in params_names:
            evaluation_point.add_argument(params_dict[name])
        evaluation_point.set_return_value(self.run_metrics_all(params_dict))
        return True

    def optimize(self, time_limit=None, reuse=False):
        import localsolver

        with localsolver.LocalSolver() as ls:
            # Create a simple model with an external function
            model = ls.model
            model_params = self.make_localsolver_parameters(model)
            f = model.create_double_external_function(
                self.evaluate_localsolver)
            func_call = model.call(f, *model_params)
            model.minimize(func_call)
            surrogate_params = f.external_context.enable_surrogate_modeling()
            model.close()

            # Setup the initial values as the default parameters
            params_dict = BenchmarkRunner.default_params()
            for name, p in zip(self.get_localsolver_param_names(), model_params):
                p.value = params_dict[name]

            if reuse:
                self.save_localsolver_solutions(surrogate_params)

            if time_limit is not None:
                ls.param.time_limit = time_limit

            # Solve
            ls.solve()

            # Get the result
            result = {}
            for name, p in zip(self.get_localsolver_param_names(), model_params):
                result[name] = p.value
            print(f"Optimization result: {result}")


parser = argparse.ArgumentParser()
parser.add_argument(
    "benchmarks", help="Benchmark instances to optimize", nargs="+", type=str)
parser.add_argument("--time-quality-tradeoff",
                    help="Tradeoff between time and quality; higher is higher quality", type=float, default=1.0)
parser.add_argument("--time-limit",
                    help="Time limit for optimization", type=int)
parser.add_argument("--reuse",
                    help="Reuse previously evaluated incumbents from the start", action="store_true")

args = parser.parse_args()

BenchmarkRunner.create_db()

runner = BenchmarkRunner(
    args.benchmarks, time_for_quality=args.time_quality_tradeoff)
runner.optimize(time_limit=args.time_limit, reuse=args.reuse)
