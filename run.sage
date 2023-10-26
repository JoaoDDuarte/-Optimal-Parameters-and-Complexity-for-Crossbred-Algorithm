import logging
import os
import warnings

# Quick and dirty way of dealing with sage imports.
os.system("sage --preparse admissibility_and_complexity.sage")
os.system("mv admissibility_and_complexity.sage.py admissibility_and_complexity.py")

os.system("sage --preparse utils.sage")
os.system("mv utils.sage.py utils.py")

os.system("sage --preparse experiments.sage")
os.system("mv experiments.sage.py experiments.py")

warnings.filterwarnings("ignore", category=DeprecationWarning)

from experiments import (
    get_fes_complexity,
    get_fes_complexity_from_1_till_n,
    get_optimal_crossbred_parameters,
    get_optimal_crossbred_params_from_1_till_n,
    get_optimal_hybrid_parameters,
    get_optimal_hybrid_params_from_1_till_n,
)
from utils import parse_args

if __name__ == "__main__":
    """Entry point for script.

    :raises Exception: n must be non-negative.
    :raises Exception: q must be larger than 1.
    """
    args = parse_args()

    algorithm = args.algorithm

    quick = args.quick
    quickest = args.quickest

    quick = True if quickest else quick

    n = args.n

    if n <= 10:
        logging.warn(
            " n <= 10, which will result in strange optimal parameters! "
            "This is a known issue. Assume D=2, d=1, k=1 to be optimal.\n"
        )
    m = args.m
    if m <= n:
        raise Exception(
            f"m={m} cannot be smaller than or equal to n={n} as the multivariate "
            "polynomial system must be overdetermined."
        )

    max_n = args.max_n
    r = args.r
    q = args.q

    min_D = args.min_D

    if n < 1:
        raise Exception("n must be larger than 0.")
    if q < 2:
        raise Exception("Finite field must be larger than 1.")
    if args.experiment:
        file_name = args.filename
        if not file_name.endswith(".csv"):
            file_name += ".csv"

        if algorithm.lower() == "crossbred":
            get_optimal_crossbred_params_from_1_till_n(
                file_name, max_n, r, q, quick=quick, quickest=quickest
            )
        elif algorithm.lower() == "hybridf5":
            get_optimal_hybrid_params_from_1_till_n(file_name, max_n, r, q)
        elif algorithm.lower() == "fes":
            get_fes_complexity_from_1_till_n(file_name, max_n, r, q)
    else:
        if algorithm.lower() == "crossbred":
            get_optimal_crossbred_parameters(n, m, q, min_D=min_D)
        elif algorithm.lower() == "hybridf5":
            get_optimal_hybrid_parameters(n, m, q)
        elif algorithm.lower() == "fes":
            get_fes_complexity(n, m, q)
