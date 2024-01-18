# Optimal Parameters and Complexity for Crossbred Algorithm


[![DOI](https://zenodo.org/badge/745179880.svg)](https://zenodo.org/doi/10.5281/zenodo.10530200)


## Description
For an overview of the Crossbred algorithm, refer to [the original paper describing it](https://hal.archives-ouvertes.fr/hal-01981516/file/2017-372.pdf).

For an overview of the Hybrid F5 algorithm, refer to [the original paper describing it](https://www.degruyter.com/document/doi/10.1515/JMC.2009.009/html).

For an overview of the FXL Algorithm, refer to the [the original paper describing it](https://arxiv.org/pdf/1112.6263v1.pdf).

For an overview of FES, refer to the [original paper descriving it](https://eprint.iacr.org/2010/313.pdf).

For the paper discussing the complexity of the Crossbred algorithm, refer to [the IACR ePrint submission](https://eprint.iacr.org/2023/1664.pdf).


This SageMath script does the following:

1. Output optimal parameters and complexity for using the Crossbred, FES, FXL or Hybrid F5 algorithm to solve a multivariate polynomial system with `n` variables, `m` equations over finite field `q >= 2`.
2. Runs an experiment which, given a relation between `m` and `n` (e.g., `m = 2n` or `m=nlog n`), for every `n` from 1 till `max_n`, calculates the optimal parameters for either the Crossbred or Hybrid-F5 algorithm. Results are displayed on screen and written to a file.

Note, if `q=2`, it is assumed that we are working in a quotient field (no squares).
Furthermore, this does not take memory considerations into account, hence optimal refers solely to computation speed.

## Prerequisites

- SageMath (tested on version 9.2, 9.8, 10.1)

This script was tested on Arch Linux with kernel 6.0.8, 6.3.2 and 6.5.8, so there may be issues with Windows and/or Mac.
However, this is unlikely.

## How to run

```
usage: run.sage.py [-h] [-n N] [-m M] [-q Q] [--min-D MIN_D] [-a {crossbred,hybridf5,fes,FXL}]
                   [--run-experiment | --no-run-experiment] [-u MAX_N] [-r R] [-o FILENAME] [--quick | --no-quick]
                   [--quickest | --no-quickest] [-d] [-v]

Script for obtaining optimal parameters for the Crossbred algorithm and run an experiment which finds the optimal
parameters for various sizes of multivariate polynomial systems. The various accepted arguments allow for different
variations of the experiment to be conducted. Note that finite field of size 2 is considered to be square-free.

options:
  -h, --help            show this help message and exit

Properties of the multivariate polynomial:
  Configure the size and field of multivariate polynomial system.

  -n N, --num-variables N
                        Number of variables. If running an experiment, If running an experiment, this will be
                        treated as the maximum value of n.
  -m M, --num-equations M
                        Number of equations. If running an experiment, this will be ignored.
  -q Q, --field Q       Finite field.
  --min-D MIN_D         Minimum value of D (ignored when running experiments as min_D is calculated per iteration).
  -a {crossbred,hybridf5,fes,FXL}, --algorithm {crossbred,hybridf5,fes,FXL}
                        Algorithm to calculate the complexity of.

Experiment:
  Configure the experiment.

  --run-experiment, --no-run-experiment
                        Whether to run the experiment.
  -u MAX_N, --up-to-n MAX_N
                        Maximum value of n.
  -r R, --relation-between-m-n R
                        Relationship between m and n, such as 2n, n+1 or nlogn.
  -o FILENAME, --filename FILENAME
                        CSV File to save results.
  --quick, --no-quick   Whether to perform an faster version of the experiment (very rarely, may produce some
                        slightly sub-optimal results). Does not apply to hybrid F5, fes and FXL.
  --quickest, --no-quickest
                        Whether to perform an fastest version of the experiment (may produce some slightly sub-
                        optimal results, takes priority over --quick). Does not apply to hybrid F5, fes and FXL.

Logging:
  Configure verbosity of logging

  -d, --debug           Print debugging statements
  -v, --verbose         Be verbose
```

### Example #1: Finding Optimal Parameters
For example, to determine the optimal parameters for solving a multivariate polynomial system with 50 variables and 60 equations over `GF(3)` using the Crossbred algorithm, run:

```python
sage run.sage -n 50 -m 60 -q 3
```

Without any verbosity or debugging options, the script will output the following:

```
Finding optimal parameters for crossbred with n=50, m=60 and q=3...

Optimal parameters found: c=142013029084110873643860, D=10, d=1, k=4
```

To speed things up, utilise the flag `--min-D` which allows you to specify a minimum value of D. 

### Example #2: Experiment

To perform the experiment from `n = 1` to `200`, `m=5n` and `q=2`, run:

```python
sage run.sage --run-experiment -u 200 -r 5n -q 2 -o experiment.csv
```

To speed things up, utilise the flag `--quick`, whereby the minimum value of `D` for `(n,m,q)` is the `D` from the optimal parameters for `(n-1,m,q)`.
Furthermore, the "not quick" algorithm will begin searching the list of sorted complexities and associated parameters from index 0. When passing the `--quick` flag,
the algorithm will take into account where the last complexity was found in the list and only begin searching from around there.  
This produces a **significant** speed up whilst sacrificing very little precision.

The reason why this works is because as `n` and `m` get bigger, there are more inadmissible parameters.

```python
sage run.sage --run-experiment -u 100 -r 2n -q 2 -o experiment.csv --quick
```

## Known issues and future improvements

For `n <= 10`, will result in strange optimal parameters. This is probably due to small missing terms in the complexity estimate. In this case, it is assumed `D=2`, `d=1`, `k=1` to be optimal.

In short, any small discrepancy in the complexity or parameter selection is most likely due to these small missing terms.

Furthermore, it may be interesting to also allow to further specify the ranges of `D`, `d`, and `k` with arguments like `min_k`, etc.

## List of Optimal Parameters

CSVs containing precomputed optimal parameters for the Crossbred, FXL and Hybrid F5 for `n=1...200` and `m=2n, n+1, nlogn` and `q=2, 3, 7` can be found in the folder named `optimal_params`.