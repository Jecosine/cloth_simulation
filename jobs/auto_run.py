import os
import sys
import itertools
import subprocess

EXEC_PATH = "./kernel_omp_papi"
# Number of nodes per dimension(n):  20, 50, 80, 110, 140, 170, 200
# Grid separation(s):                0.1, 0.5, 1.0, 1.5, 2.0
# Mass of node(m):                   0.1, 0.5, 1.0, 1.5, 2.0, 5.0, 10.0
# Force constant(f)                  1.0, 5.0, 10.0, 20.0, 50.0, 100.0
# Node Interaction Level (delta): 2, 3, 4
# Gravity(g):                        0.981000
# Radius of Ball(b):                 1.0, 3.0, 5.0, 10.0
# Offset of falling cloth(o):        3.0
# Timestep(t):                       0.01, 0.02, 0.05, 0.1, 0.5, 1
# num iterations(i):                 400, 500
ns = [200]
ss = [1.0]
ms = [1.0]
fs = [10.0]
ds = [i for i in range(1, 11)]
gs = [0.981]
bs = [3.0]
oos = [3.0]
ts = [0.05]
iis = [400]
ps = [16]

# param_list = [ns, ss, ms, fs, ds, gs, bs, oos, ts, iis]
# param_names = ["n", "s", "m", "f", "d", "g", "b", "o", "t", "i"]
# # omp version
param_list = [ns, ss, ms, fs, ds, gs, bs, oos, ts, iis, ps]
param_names = ["n", "s", "m", "f", "d", "g", "b", "o", "t", "i", "p"]
assert len(param_list) == len(param_names)
def run_with_params(param_names, values, f):
    args = sum([[f"-{p}",str(v)] for p, v in zip(param_names, values)], [])
    cmd = [EXEC_PATH, *args]
    print(f"Running with cmd: {cmd}")
    subprocess.Popen(cmd, stdout=f).wait()

def check_executable(path: str):
    if not os.path.exists(path):
      raise RuntimeError(f"executable {path} not exists")
    return True


if __name__ == "__main__":
  print("------ Running executable: origin ------")
  check_executable(EXEC_PATH)
  all_combs = list(itertools.product(*param_list))
  print(f"testing {len(all_combs)} combinations...")
  counter = 0
  with open("omp_p.log", "a") as f:
    for comb in all_combs:
        print(f"Runing {comb} ...")
        for i in range(10):
          run_with_params(param_names, comb, f)
  print("All done")


