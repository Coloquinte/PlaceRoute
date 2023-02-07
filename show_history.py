
import pickle
import numpy as np
import matplotlib.pyplot as plt
import sys


history_file = "history.pkl"
if len(sys.argv) >= 2:
    history_file = sys.argv[1]

with open(history_file, 'rb') as f:
    history = pickle.load(f)

def get_mean(history):
    main_keys = ["benchmark", "global.seed", "detailed.seed"]
    metrics = ["metrics.time_total", "metrics.time_global", "metrics.time_detailed", "metrics.hpwl"]
    benchmarks = []
    for h in history:
        key_vals = [h[k] for k in main_keys]
        if key_vals not in benchmarks:
            benchmarks.append(key_vals)
    unique_hist = []
    for h in history:
        uh = dict(h)
        for k in main_keys:
            del uh[k]
        for k in metrics:
            del uh[k]
        if uh not in unique_hist:
            unique_hist.append(uh)
    ret = []
    for uh in unique_hist:
        metric_vals = {m : [] for m in metrics}
        missing_data = False
        for key_vals in benchmarks:
            h = dict(uh)
            for k, v in zip(main_keys, key_vals):
                h[k] = v
            datapoint = None
            for o in history:
                found = True
                for k, v in h.items():
                    if o[k] != v:
                        found = False
                if found:
                    datapoint = o
                    break
            if datapoint is None:
                print(f"Missing benchmark for {uh}: {key_vals}")
                missing_data = True
            else:
                for m in metrics:
                    metric_vals[m].append(datapoint[m])
        if missing_data:
            continue
        val = dict(uh)
        for m in metrics:
            mean = np.exp(np.log(metric_vals[m]).mean())
            val[m] = mean
        ret.append(val)
    return ret

def sorted_history(history):
    return sorted(history, key=lambda d: d["metrics.hpwl"])

def get_pareto(history):
    tuples = []
    for h in history:
        tuples.append((h["metrics.time_total"], h["metrics.hpwl"]))
    non_dominated = []
    for h1 in history:
        ok = True
        for h2 in history:
            if h2["metrics.time_total"] < h1["metrics.time_total"] and h2["metrics.hpwl"] < h1["metrics.hpwl"]:
                ok = False
        if ok:
            non_dominated.append(h1)
    return sorted(non_dominated, key=lambda d: d["metrics.hpwl"])

def show_metrics(pareto, name):
    metrics = [h[name] for h in pareto]
    if (np.array(metrics) == metrics[0]).all():
        return
    plt.title(f"Pareto front: {name}")
    plt.plot(metrics)
    plt.show()

def show_cross(history, pareto, m1, m2):
    if m1 == m2:
        return
    v1 = [h[m1] for h in pareto]
    v2 = [h[m1] for h in pareto]
    if (np.array(v1) == np.array(v1)[0]).all():
        return
    if (np.array(v2) == np.array(v2)[0]).all():
        return
    plt.title(f"Comparison {m1} {m2}")
    v1 = [h[m1] for h in history]
    v2 = [h[m2] for h in history]
    plt.scatter(v1, v2)
    v1 = [h[m1] for h in pareto]
    v2 = [h[m2] for h in pareto]
    plt.scatter(v1, v2, color="red")
    plt.show()

def show_pareto(history, pareto):
    times = [h["metrics.time_total"] for h in history]
    qualities = [h["metrics.hpwl"] for h in history]
    plt.scatter(times, qualities)
    times = [h["metrics.time_total"] for h in pareto]
    qualities = [h["metrics.hpwl"] for h in pareto]
    plt.title(f"Pareto front")
    plt.scatter(times, qualities, color="red")
    plt.show()

def print_params(pareto):
    names = list(pareto[0].keys())
    used_names = []
    for name in names:
        metrics = [h[name] for h in pareto]
        if (np.array(metrics) == metrics[0]).all():
            continue
        used_names.append(name)
    names = used_names
    used_params = [23, 22, 18, 16, 14, 12, 9, 5, 1]
    for name in names:
        l = [pareto[i][name] for i in used_params]
        print(f"\t{name}: {l}")

history = get_mean(history)
pareto = get_pareto(history)

show_pareto(history, pareto)
names = list(pareto[0].keys())
for name in names:
    show_metrics(pareto, name)
for name in names:
    show_cross(history, pareto, "metrics.hpwl", name)

