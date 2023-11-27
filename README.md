# Simulation

## Files Descriptions

- `data/`: This directory contains the csv data recording the performance of program. The format of file please see "Data format" section.(The data will not be public to github)
  - `main.csv`: performance of `main`, vary $n$ from $20$ to $300$ with step $20$
  - `opt.csv`: performance of `opt`, vary $n$ from $20$ to $300$ with step $20$
  - `vect_omp.csv`: performance of `vect_omp`, vary $n$ from $20$ to $300$ with step $20$
  - `sse.csv`: performance of `sse`, vary $n$ from $20$ to $300$ with step $20$
  - `omp_with_threads.csv`: performance of `omp`, $n=1000$, $p$ varies from $1$ to $24$
  - `omp_with_threads_48.csv`: performance of `omp`, $n=1000$, $p$ varies from $1$ to $48$, run on compute node.
  - `omp_1000_2.csv`, `omp_500_4`, `omp_250_8`: performance of `omp`, $p=24$, $(n, d) = \{(1000, 2), (500, 4), (250, 8)\}$
  - `omp_chunk_default.csv`, `omp_chunk_1.csv`, `omp_chunk_10.csv`: performance of `omp`, $p$ varies from $1$ to $24$ for each chunksize$=\{\text{default}, 1, 10\}$
  - `comp6464_p.csv`, `omp_p.csv`: performance of `omop_comp6464` and `omp`, $d$ varies from $1$ to $10$ with fixed $n=200, p=16$
- `plots/`: This directory contains plots generated based on csv from `data/`
  - `plot.ipynb`: Jupyter notebook that produce plots for `writeup.pdf` and other figures in this directory.
  - `roofline_omp.html`: roofline plot
- `jobs/`: contains the batch job submission to compute node
  - `batch_job`: job description
  - `auto_run.py`: run specified parameter combinations and output to a file
- `*.[cpp|h]`: See `CMakeLists.txt` for dependencies between executables and code files
- `CMakeLists.txt`


## Usage

On gadi, fisrtly, load all necessary modules

```bash
module load papi intel-compiler cmake/3.18.2 python3-as-python
```

Configure cmake, ensure you are at the project root folder

```bash
mkdir build && cd build && rm -rf ./* && cmake ..
```

From the output log, confirm you can see `PAPI FOUND` which is needed by some executables measuring the performance.

You can make all executables or one of them (Now you are in `build/` directory):
```bash
# build all exe
make -j24
# OR just build one of them
make <target_name>
```

The available targets are listed below:

- `kernel_main`: the main version
- `kernel_opt`: the opt version
- `kernel_vect_omp`: the version vectorized by OpenMP
- `kernel_sse`: the version vectorized by SSE (AVX256)
- `kernel_omp`: the omp version
- `kernel_omp_comp6464`: the modified block-wised omp version required by COMP6464
- `kernel_main_papi`: the papi version output metrics for main
- `kernel_opt_papi`: the papi version output metrics for opt
- `kernel_vect_omp_papi`: the papi version output metrics for vect_omp
- `kernel_sse_papi`: the papi version output metrics for sse
- `kernel_omp_papi`: the papi version output metrics for omp
- `kernel_omp_comp6464_papi`: the papi version output metrics for omp_comp6464

## Data Format

The output of papi executable are:

```python
[
  "n", 
  "sep", 
  "mass", 
  "fcon", 
  "delta", 
  "g", 
  "rball", 
  "offset", 
  "dt", 
  "maxiter", 
  "time",
  "PAPI_TOT_INS", 
  "PAPI_TOT_CYC",
  "PAPI_L1_DCM",
  "PAPI_L3_TCA",
  "PAPI_L3_TCM",
  "PAPI_L2_DCM",
  "PAPI_L2_DCA",
  "PAPI_BR_MSP",
  "PAPI_BR_CN"
]
```

The `kernel_omp_papi` provide `double` type "time", while the others provide `long`

## Optimization Insight

- [Serial optimization version](./docs/serial.md)

- [SSE(AVX256) version](./doc/sse.md)

- [Roofline analysis](./doc/roofline.md)