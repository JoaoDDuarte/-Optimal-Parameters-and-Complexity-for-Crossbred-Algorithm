import argparse
import logging
import os
import re

from sage.functions.log import log
from sage.functions.other import ceil
from sage.rings.all import ZZ
from sage.rings.polynomial.polydict import ETuple

os.system("sage --preparse admissibility_and_complexity.sage")
os.system("mv admissibility_and_complexity.sage.py admissibility_and_complexity.py")

from admissibility_and_complexity import crossbred_complexity, degree_of_regularity


def is_admissible_crossbred(admis, D, d, k):
    """
    Determines whether parameters D, d and k are admissible.

    :param admis: Admissibility generating series.
    :param D: Degree of Macaulay Matrix,
    :param d: The desired degree of the system of multivariate polynomials after
        we specialise the last n - k variables.
    :param k: Number of variables we want our specialised system of multivariate
        polynomials to have.
    :return: True if parameters are admissible.
    """
    d_tuple = ETuple([D, d])
    logging.debug(f"\t\tCoefficient of {admis[d_tuple]} for D={D}, d={d}, k={k}")
    return admis[d_tuple] >= 0


def m_to_n(relation, n):
    """
    Converts a string indicating the relationship between m and n
    to an actual number.

    :param relation: Relationship between m and n, e.g. 2n.
    :param n: Number of equations.
    :raises Exception: Invalid or malformed relation between m and n.
    :return: The evaluation of the relationship, e.g., if n = 3 and m = 3n,
        then return 6.
    """
    relation = relation.replace(" ", "")

    r_multiplication = re.compile("[0-9]+n\Z")
    r_addition = re.compile("n\+[0-9]+\Z")

    if r_multiplication.match(relation) is not None:
        return int(relation.replace("n", "")) * n
    elif r_addition.match(relation):
        return n + int(relation.replace("n", "").replace("+", ""))
    elif relation == "nlogn":
        return n * ZZ(ceil(log(n, 10)))
    else:
        raise Exception(
            f"Relation {relation} between n and m not supported! Use (e.g.) 2n, n+1 or nlogn."
        )


def get_complexity_list_and_dictionary(n, m, q, min_k, min_d, min_D):
    """
    Returns complexity list and dictionary for the crossbred complexity analysis.

    :param n: Number of variables.
    :param m: Number of equations.
    :param min_k: Minimum value of k, defaults to -1.
    :param min_d: Minimum value of d, defaults to 1.
    :param min_D: Minimum value of D, defaults to -1.
    :return: Complexity list and dictionary.
    """

    deg_n_m = degree_of_regularity(n, m, q)

    complexity_dictionary = {}
    precomputed_overall_min_D = degree_of_regularity(min_k, m, q)
    if min_D != -1 and min_D < precomputed_overall_min_D:
        precomputed_overall_min_D = min_D
    logging.info(
        f"\tCalculating all possible admissible and inadmissible parameters with minimum D = {precomputed_overall_min_D}"
        f" and max D={deg_n_m}..."
    )

    for k in range(min_k, n + 1):
        deg_k_m = degree_of_regularity(k, m, q)
        minimum_D = (
            deg_k_m if (min_D == -1) or (min_D != -1 and deg_k_m > min_D) else min_D
        )
        for D in range(minimum_D, deg_n_m + 1):
            for d in range(min_d, D):
                complexity_dictionary[crossbred_complexity(D, d, n, k, q)] = (D, d, k)

    complexity_list = list(complexity_dictionary.keys())
    complexity_list.sort()
    return complexity_dictionary, complexity_list


def get_starting_position(q, start_searching_at_index, complexity_list_len):
    """
    Get starting position to start searching a list of crossbred complexities for the optimal params/complexity.

    :param q: Finite field in which the system of multivariate polynomials are defined.
    :param start_searching_at_index: Where to begin searching for the minimum in the complexity list.
    :param complexity_list_len: Length of the complexity list.
    :return: Where to begin searching for the minimum in the complexity list.
    """
    start_searching_at_index = (
        start_searching_at_index - 10 if q == 2 else start_searching_at_index - 20
    )
    start_searching_at_index = (
        0 if start_searching_at_index < 0 else start_searching_at_index
    )

    if complexity_list_len < start_searching_at_index:
        logging.debug(
            f"\tPosition at where to start searching = {start_searching_at_index} is larger than {complexity_list_len} - "
            "setting start position to 0..."
        )
        start_searching_at_index = 0
    return start_searching_at_index


def parse_args():
    """
    Parse arguments.

    :return: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Script for obtaining optimal parameters for the Crossbred algorithm "
        "and run an experiment which finds the optimal parameters for various sizes of "
        "multivariate polynomial systems. The various accepted arguments allow for different"
        " variations of the experiment to be conducted. Note that finite field of size 2 is "
        "considered to be square-free."
    )

    group_optimal_params = parser.add_argument_group(
        "Properties of the multivariate polynomial",
        "Configure the size and field of multivariate polynomial system.",
    )
    group_optimal_params.add_argument(
        "-n",
        "--num-variables",
        help="Number of variables. If running an experiment, If running an experiment, "
        "this will be treated as the maximum value of n.",
        type=int,
        dest="n",
        default=30,
    )
    group_optimal_params.add_argument(
        "-m",
        "--num-equations",
        help="Number of equations. If running an experiment, this will be ignored.",
        type=int,
        dest="m",
        default=60,
    )
    group_optimal_params.add_argument(
        "-q",
        "--field",
        help="Finite field.",
        type=int,
        dest="q",
        default=2,
    )
    group_optimal_params.add_argument(
        "--min-D",
        help="Minimum value of D (ignored when running experiments as min_D is calculated per iteration).",
        type=int,
        dest="min_D",
        default=-1,
    )
    group_optimal_params.add_argument(
        "-a",
        "--algorithm",
        help="Algorithm to calculate the complexity of.",
        type=str,
        choices=["crossbred", "hybridf5", "fes", "FXL"],
        dest="algorithm",
        default="crossbred",
    )
    group_experiment = parser.add_argument_group(
        "Experiment",
        "Configure the experiment.",
    )
    group_experiment.add_argument(
        "--run-experiment",
        action=argparse.BooleanOptionalAction,
        dest="experiment",
        default=False,
        help="Whether to run the experiment.",
    )
    group_experiment.add_argument(
        "-u",
        "--up-to-n",
        help="Maximum value of n.",
        type=int,
        dest="max_n",
        default=100,
    )
    group_experiment.add_argument(
        "-r",
        "--relation-between-m-n",
        help="Relationship between m and n, such as 2n, n+1 or nlogn.",
        type=str,
        dest="r",
        default="2n",
    )
    group_experiment.add_argument(
        "-o",
        "--filename",
        dest="filename",
        help="CSV File to save results.",
        type=str,
        default="experiment.csv",
    )
    group_experiment.add_argument(
        "--quick",
        action=argparse.BooleanOptionalAction,
        dest="quick",
        default=False,
        help="Whether to perform an faster version of the experiment"
        " (very rarely, may produce some slightly sub-optimal results)."
        " Does not apply to hybrid F5, fes and FXL.",
    )
    group_experiment.add_argument(
        "--quickest",
        action=argparse.BooleanOptionalAction,
        dest="quickest",
        default=False,
        help="Whether to perform an fastest version of the experiment"
        " (may produce some slightly sub-optimal results, takes priority"
        " over --quick). Does not apply to hybrid F5, fes and FXL.",
    )
    group_logging = parser.add_argument_group(
        "Logging", "Configure verbosity of logging"
    )
    group_logging.add_argument(
        "-d",
        "--debug",
        help="Print debugging statements",
        action="store_const",
        dest="loglevel",
        const=logging.DEBUG,
        default=logging.WARNING,
    )
    group_logging.add_argument(
        "-v",
        "--verbose",
        help="Be verbose",
        action="store_const",
        dest="loglevel",
        const=logging.INFO,
    )

    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel)

    return args
