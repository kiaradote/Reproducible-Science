"""Stress run of the notebook in a fresh process with a memory guard (stops before the machine starts swapping)."""
import os, sys, time, threading, json
sys.path.insert(0, os.getcwd())
import psutil, run_nb
dataset, limit_mb = sys.argv[1], int(sys.argv[2])
state = {"after_cell": None, "peak": 0}
vm = psutil.virtual_memory()
print(f"start: total {vm.total/2**30:.1f} GiB, available {vm.available/2**30:.2f} GiB, guard {limit_mb} MB RSS", flush=True)
def watch():
    p = psutil.Process()
    while True:
        rss = p.memory_info().rss / 2**20
        state["peak"] = max(state["peak"], rss)
        if rss > limit_mb:
            print(f"MEMORY GUARD: RSS {rss:.0f} MB exceeded {limit_mb} MB after cell {state['after_cell']}; stopping to protect the machine", flush=True)
            os._exit(3)
        time.sleep(0.5)
threading.Thread(target=watch, daemon=True).start()
t0 = time.time()
log, ns = run_nb.run_notebook("scrna_checkpoints.ipynb", env={"NB_DATASET": dataset},
                              after_cell=lambda i, n: (state.update(after_cell=i), print(f"[{time.time()-t0:7.1f}s] cell {i} done, RSS now {psutil.Process().memory_info().rss/2**20:.0f} MB, peak {state['peak']:.0f} MB", flush=True)))
for e in log:
    print(f"cell {e['cell']}: {'ok' if e['ok'] else 'ERROR ' + str(e['error'])[:300]}")
print("FINAL PEAK RSS MB", round(state["peak"]), "seconds", round(time.time() - t0))
print(ns["wf"].status().to_string())
