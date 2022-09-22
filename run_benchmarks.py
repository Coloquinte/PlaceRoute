import math
import os
import sqlite3
import time

import coloquinte

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
        circuit = coloquinte.Circuit.read_ispd(os.path.join(self._prefix, self.benchmark))
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
    def __init__(self):
        self.con = sqlite3.connect("benchmarks.db")

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
        unique_columns = list(k for k in data.keys() if k.startswith("global") or k.startswith("detailed"))
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
        cols = BenchmarkRunner.unique_columns(data)
        column_names = []
        columns = []
        values = []
        for k, v in data.items():
            column_names.append(k)
            columns.append(k + "=?")
            values.append(v)
        columns = " and ".join(columns)
        sql = f"SELECT * FROM ColoquinteBenchmarks WHERE {columns}"
        c = self.con.cursor()
        c.execute(sql, values)
        ret = c.fetchall()
        if len(ret) == 0:
            return None
        elif len(ret) == 1:
            return {k: v for k, v in zip(column_names, ret[0])}
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
        c = self.con.cursor()
        c.execute(sql, values)
        self.con.commit()
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
        return data

    def run_metrics(self, params_dict, time_for_quality):
        """
        Return a metrics of the run, with a tradeoff between time and quality
        """
        assert time_for_quality > 0
        m = self.run(params)
        t = m["time_total"]
        q = m["hpwl"]
        return math.log(q) - time_for_quality * math.log(t)

    def run_metrics_all(self, params_dict, time_for_quality, benchmarks=None, prefix="benchmarks/ISPD06/"):
        """
        Return the geometric mean of the metrics across the benchmarks, with the given time/quality tradeoff
        """
        if benchmarks is None:
            benchmarks = os.listdir(prefix)
        metrics = []
        for benchmark in sorted(benchmarks):
            params = dict(params_dict)
            params["benchmark"] = benchmark
            params["prefix"] = prefix
            metrics.append(self.run_metrics(params))
        return np.exp(np.mean(np.log(metrics)))

    def optimize(self, benchmarks=None):
        import localsolver

        with localsolver.LocalSolver() as ls:
            model = ls.model
            max_nb_steps = ls.int(20, 60)
            gap_tolerance = ls.float(0.01, 0.2)


BenchmarkRunner.create_db()

runner = BenchmarkRunner()
runner.optimize(["adaptec1"])
