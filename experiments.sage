import logging
import os

os.system("sage --preparse admissibility_and_complexity.sage")
os.system("mv admissibility_and_complexity.sage.py admissibility_and_complexity.py")

os.system("sage --preparse utils.sage")
os.system("mv utils.sage.py utils.py")

from sage.functions.log import log
from sage.functions.other import floor
from sage.rings.all import ZZ

from admissibility_and_complexity import (
    crossbred_admissibility_series,
    hybrid_f5_complexity,
)
from utils import (
    get_complexity_list_and_dictionary,
    get_starting_position,
    is_admissible_crossbred,
    m_to_n,
)


def get_optimal_hybrid_parameters(n, m, q):
    """Calculate the optimal parameters for the Hybrid F5 algorithm.

    :param n: Number of variables.
    :param m: Number of equations.
    :param q: Finite field in which the system of multivariate polynomials are defined.
    :return: Complexity and optimal parameters (in this case, value of k).
    """
    print(f"Finding optimal parameters for hybrid f5 with n={n}, m={m} and q={q}...\n")

    complexity = {}
    try:
        for k in range(1, n):
            c = hybrid_f5_complexity(n, m, k, q)
            complexity[c] = k
        complexity_list = list(complexity.keys())
        min_c = min(complexity_list)
        optimal_k = complexity[min_c]

        print(f"Optimal parameters found: c={min_c} and k={optimal_k}\n")
        return min_c, optimal_k
    except (TypeError, ValueError) as e:
        print(e)
        return q, 1


def get_fes_complexity(n, m, q):
    """Calculate the complexity for FES.

    :param n: Number of variables.
    :param m: Number of equations.
    :param q: Finite field in which the system of multivariate polynomials are defined.
    :return: Complexity and optimal parameters (in this case, value of k).
    """
    print(f"Finding complexity for FES with n={n}, m={m} and q={q}...\n")
    c = ZZ(floor(log(n, q))) * q**n
    print(f"Complexity for FES calculated: c = {c}\n")
    return c


def get_optimal_crossbred_parameters(
    n, m, q, min_D=-1, min_d=1, min_k=1, start_searching_at_index=0
):
    """Calculate the optimal parameters for the Crossbred algorithm.
    The optimisation comes from not checking values of D below a certain
    number. These can be circumvented by setting min_d = 1, min_k = 1,
    min_D = -1 and max_D = -1.

    :param n: Number of variables.
    :param m: Number of equations.
    :param q: Finite field in which the system of multivariate polynomials are defined.
    :param min_D: Minimum value of D, defaults to -1.
    :param min_d: Minimum value of d, defaults to 1.
    :param min_k: Minimum value of k, defaults to -1.
    :param start_searching_at_index: Where to begin searching for the minimum in the complexity list.
    :return: Complexity and optimal parameters.
    """

    print(f"Finding optimal parameters for crossbred with n={n}, m={m} and q={q}...\n")

    min_d = min_d if min_d >= 1 else 1
    min_k = min_k if min_k >= 1 else 1

    complexity_dictionary, complexity_list = get_complexity_list_and_dictionary(
        n, m, q, min_k, min_d, min_D
    )
    start_searching_at_index = get_starting_position(
        q, start_searching_at_index, len(complexity_dictionary)
    )
    complexity_list = complexity_list[start_searching_at_index:]

    logging.info("\tChecking admissibility of all possible parameters...")

    where_in_list = start_searching_at_index - 1
    for c in complexity_list:
        where_in_list += 1
        D, d, k = complexity_dictionary[c]
        logging.debug("")
        logging.debug(
            f"\t\tAdmissibility series is being calculated for n={n}, m={m}, q={q}"
            f" c={c} D={D}, d={d}, k={k}..."
        )
        admiss_series = crossbred_admissibility_series(n, m, k, q, D, d).dict()
        logging.debug("\t\tDone.")
        try:
            if n <= 10:
                logging.debug("\t\t=> Admissible.\n")
                print(f"Optimal parameters found: c={c}, D=2, d=1, k=1\n")
                return c, 2, 1, 1, 0

            if is_admissible_crossbred(admiss_series, D, d, k):
                c_D, c_d, c_k = complexity_dictionary[c]
                logging.debug(f"\t\t=> Admissible result at index {where_in_list}.\n")
                print(f"Optimal parameters found: c={c}, D={c_D}, d={c_d}, k={c_k}\n")

                return c, c_D, c_d, c_k, where_in_list
            else:
                logging.debug(
                    "\t\t=> Inadmissible, trying parameters that produce the next best complexity..."
                )
        except Exception as e:
            logging.error(
                f" No optimal parameters found for n={n}, m={m} and q={q}, produced error {e}\n"
            )


def get_optimal_hybrid_params_from_1_till_n(name, max_n, relation, q):
    """Given a relation between m and n (e.g., m = 2n), for every n from 1 till max_n,
    calculate the optimal parameters for hybrid F5. Results are displayed
    on screen and written to a file. If a file is present with the same name, it will
    append the results.

    :param name: Name of the file to write the results.
    :param max_n: Maximum value of n.
    :param q: Finite field in which the system of multivariate polynomials are defined.
    :param relation: Relationship between m and n, e.g. 2n.
    """
    exists = os.path.exists(name)
    empty = False

    start = 1

    if exists:
        if os.stat(name).st_size == 0:
            empty = True
        else:
            start = sum(1 for _ in open(name))
            print("Resuming experiment...")

    if relation.replace(" ", "") == "nlogn":
        start += 10
    logging.debug(f"\tStarting from n={start}\n")

    with open(name, "a") as output_file:
        if empty or not exists:
            output_file.write("n,m,c,k\n")
        for n in range(start, max_n + 1):
            m = m_to_n(relation, n)
            c, k = get_optimal_hybrid_parameters(n, m, q)
            output_file.write(f"{n},{m},{c},{k}\n")
            logging.debug(f"\tWrote values n={n},m={m},c={c},k={k} to file.\n")


def get_fes_complexity_from_1_till_n(name, max_n, relation, q):
    """Given a relation between m and n (e.g., m = 2n), for every n from 1 till max_n,
    calculate the complexity for FES. Results are displayed
    on screen and written to a file. If a file is present with the same name, it will
    append the results.

    :param name: Name of the file to write the results.
    :param max_n: Maximum value of n.
    :param q: Finite field in which the system of multivariate polynomials are defined.
    :param relation: Relationship between m and n, e.g. 2n.
    """
    exists = os.path.exists(name)
    empty = False

    start = 1

    if exists:
        if os.stat(name).st_size == 0:
            empty = True
        else:
            start = sum(1 for _ in open(name))
            print("Resuming experiment...")

    if relation.replace(" ", "") == "nlogn":
        start += 10
    logging.debug(f"\tStarting from n={start}\n")

    with open(name, "a") as output_file:
        if empty or not exists:
            output_file.write("n,m,c\n")
        for n in range(start, max_n + 1):
            m = m_to_n(relation, n)
            c = get_fes_complexity(n, m, q)
            output_file.write(f"{n},{m},{c}\n")
            logging.debug(f"\tWrote values n={n},m={m},c={c} to file.\n")


def get_optimal_crossbred_params_from_1_till_n(
    name, max_n, relation, q, quick=False, quickest=False
):
    """Given a relation between m and n (e.g., m = 2n), for every n from 1 till max_n,
    calculate the optimal parameters for crossbred. Results are displayed
    on screen and written to a file. If a file is present with the same name, it will
    append the results.

    TODO: separate into smaller functons, this is a bit monolithic

    :param name: Name of the file to write the results.
    :param max_n: Maximum value of n.
    :param relation: Relationship between m and n, e.g. 2n.
    :param q: Finite field in which the system of multivariate polynomials are defined.
    :param quick: Whether to run the quick version of the get_optimal_crossbred_parameters() algorithm.
    :param quickest: Whether to run the quickest version of the get_optimal_crossbred_parameters() algorithm.
    """
    exists = os.path.exists(name)
    empty = False

    start = 1
    min_D = -1
    min_k = 1
    min_d = 1
    last_position_in_complexity_list = 0

    if exists:
        if os.stat(name).st_size == 0:
            empty = True
        else:
            start = sum(1 for _ in open(name))
            print("Resuming experiment...")

        if quick and not empty:
            with open(name, "r") as output_file:
                last_line = output_file.readlines()[-1]

                last_position_in_complexity_list = (
                    last_position_in_complexity_list
                    if start == 1
                    else int(last_line.split(",")[6])
                )

                if quickest:
                    min_D = min_D if start == 1 else int(last_line.split(",")[3])
                    min_d = min_d if start == 1 else int(last_line.split(",")[4]) - 3
                    min_k = min_k if start == 1 else int(last_line.split(",")[5]) - 20
                    logging.debug(f"\tMinimum value of D is set to {min_D} from file.")
                    logging.debug(
                        f"\tMinimum value of d is set to {min_d} from file (if value is negative, it will default to 1)."
                    )
                    logging.debug(
                        f"\tMinimum value of k is set to {min_k} from file (if value is negative, it will default to 1)."
                    )

                logging.debug(
                    f"\tWill only search complexity list starting at index {last_position_in_complexity_list} from file."
                )
    if relation.replace(" ", "") == "nlogn":
        start += 10
    logging.debug(f"\tStarting from n={start}\n")

    with open(name, "a") as output_file:
        if empty or not exists:
            output_file.write("n,m,c,D,d,k,index_in_list\n")
        for n in range(start, max_n + 1):
            m = m_to_n(relation, n)
            try:
                c, D, d, k, l = get_optimal_crossbred_parameters(
                    n,
                    m,
                    q,
                    min_D=min_D,
                    min_d=min_d,
                    min_k=min_k,
                    start_searching_at_index=last_position_in_complexity_list,
                )
            except TypeError as _:
                logging.info(
                    "Could not populate the complexity list, rerunning current iteration without the --quick flag..."
                )
                c, D, d, k, l = get_optimal_crossbred_parameters(n, m, q)
                min_D = D
            if quick:
                last_position_in_complexity_list = l
                logging.debug(
                    f"\tMinimum value of D for next iteration is set to {min_D}."
                )
                if quickest:
                    min_D = D
                    min_d = d - 3
                    min_k = k - 20
                    logging.debug(
                        f"\tMinimum value of d for next iteration is set to {min_d} (if value is negative, it will default to 1)."
                    )
                    logging.debug(
                        f"\tMinimum value of k for next iteration is set to {min_k} (if value is negative, it will default to 1)."
                    )
                logging.debug(
                    f"\tNew index to start search for next iteration is set to {last_position_in_complexity_list}."
                )

            output_file.write(f"{n},{m},{c},{D},{d},{k},{l}\n")
            logging.debug(
                f"\tWrote values n={n},m={m},c={c},D={D},d={d},k={k},index_in_list={l} to file.\n"
            )
