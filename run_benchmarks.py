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
    def columns(data=None):
        if data is None:
            data = BenchmarkRun("").to_dict()
        columns = list(data.keys())
        return BenchmarkRunner.metric_columns(data) + columns

    @staticmethod
    def metric_columns(data=None):
        return ["time_total", "time_global", "time_detailed", "hpwl"]

    @staticmethod
    def unique_columns(data=None):
        if data is None:
            data = BenchmarkRun("").to_dict()
        unique_columns = list(k for k in data.keys() if k.startswith(
            "global") or k.startswith("detailed"))
        return ["benchmark", "prefix"] + unique_columns

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

    def get_data(self, data):
        key_columns = BenchmarkRunner.unique_columns(data)
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
        return self.time_for_quality * math.log(q) + math.log(t)

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
            "global_penalty_cutoff_distance",
            "global_approximation_distance",
            "detailed_nb_passes",
            "detailed_local_search_nb_neighbours",
        ]

    def make_localsolver_parameters(self, model):
        return [
            model.int(20, 60),
            model.float(0.01, 0.2),
            model.float(2.0, 50.0),
            model.float(0.1, 10.0),
            model.int(1, 3),
            model.int(1, 6),
        ]

    def optimize(self):
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

            # Add limits
            surrogate_params.evaluation_limit = 1000
            ls.param.time_limit = 600

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

args = parser.parse_args()

BenchmarkRunner.create_db()

runner = BenchmarkRunner(
    args.benchmarks, time_for_quality=args.time_quality_tradeoff)
runner.optimize()
