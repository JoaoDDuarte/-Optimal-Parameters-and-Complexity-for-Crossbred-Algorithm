from sage.arith.all import binomial
from sage.rings.all import QQ, ZZ
from sage.rings.power_series_ring import PowerSeriesRing


def degree_of_regularity(n, m, q):
    """Calculate the degree of regularity.

    :param n: Number of variables.
    :param m: Number of equations.
    :param q: Finite field in which the system of multivariate polynomials are defined.
    :raises Exception: Number of variables must be smaller than the number of equations.
    :return: Degree of regularity.
    """

    degs = [q for _ in range(0, m)]
    if m <= n:
        raise ValueError(
            f"This function requires an overdefined system of polynomials (n={n}, m={m})."
        )

    from sage.misc.misc_c import prod
    from sage.rings.power_series_ring import PowerSeriesRing
    from sage.rings.rational_field import QQ, ZZ

    z = PowerSeriesRing(QQ, "z", default_prec=sum(degs)).gen()
    if q == 2:
        s = (1 + z) ** n / prod([1 + z**d for d in degs])
    else:
        s = prod([1 - z**d for d in degs]) / (1 - z) ** n

    for dreg in range(sum(degs)):
        if s[dreg] <= 0:
            return ZZ(dreg)
    raise ValueError("BUG: Could not compute the degree of semi-regularity")


def crossbred_admissibility_series(n, m, k, q, D, d):
    """Returns expanded (up to precision D + d + k) of the admissibility
    generating series. It is assumed that the degree of the polynomial is q.

    :param n: Number of variables.
    :param m: Number of equations.
    :param k: Number of variables we want our specialised system of multivariate
        polynomials to have.
    :param q: Finite field in which the system of multivariate polynomials are defined.
    :param D: Degree of Macaulay Matrix,
    :param d: The desired degree of the system of multivariate polynomials after we
        specialise the last n - k variables.
    :return: Expanded admissibility generating series.
    """
    if n >= m:
        raise ValueError(f"n >= m ({n} >= {m}), system must be overdetermined")
    R = PowerSeriesRing(QQ, "X, Y", default_prec=D + d + 1)
    (X, Y) = R.gens()

    if q == 2:
        elem0 = ((1 + X) ** (n - k)) / ((1 - X) * (1 - Y))

        elem1 = ((1 + X * Y) ** k) / ((1 + X**2 * Y**2) ** m)
        elem2 = ((1 + X) ** k) / ((1 + X**2) ** m)

        S_dk = elem0 * (elem1 - elem2)

        resultant_series = S_dk - (
            ((1 + Y) ** k) / ((1 - X) * (1 - Y) * (1 + Y**2) ** m)
        )
    else:
        elem0 = (1) / ((1 - X) * (1 - Y) * (1 - X) ** (n - k))

        elem1 = ((1 - X**q * Y**q) ** m) / ((1 - X * Y) ** k)
        elem2 = ((1 - X**q) ** m) / ((1 - X) ** k)

        S_dk = elem0 * (elem1 - elem2)

        resultant_series = S_dk - (
            (1 - Y**q) ** m / ((1 - Y) ** k) * (1 - X) * (1 - Y)
        )

    return resultant_series


def crossbred_complexity(D, d, n, k, q):
    """Calculates complexity of the Crossbred algorithm.

    :param D: Degree of Macaulay Matrix,
    :param d: The desired degree of the system of multivariate polynomials after we
        specialise the last n - k variables.
    :param n: Number of variables.
    :param k: Number of variables we want our specialised system of multivariate
        polynomials to have.
    :param q: Finite field in which the system of multivariate polynomials are defined.
    :return: Complexity of the Crossbred algorithm.
    """
    f_k = 0
    for d_k in range(d + 1, D + 1):
        for d_s in range(0, D - d_k + 1):
            if q != 2:
                f_k = f_k + (
                    (binomial(k + d_k - 1, d_k) * binomial(n - k + d_s - 1, d_s))
                )
            else:
                f_k = f_k + (binomial(k, d_k) * binomial(n - k, d_s))
    if q != 2:
        linearisation = binomial(k + d - 1, d) ** 2
    else:
        linearisation = sum([binomial(k, i) for i in range(0, d + 1)]) ** 2
    return f_k**2 + q ** (n - k) * linearisation


def hybrid_f5_complexity(n, m, k, q):
    """
    Returns complexity of the Hybrid F5 complexity.

    :param n: Number of variables.
    :param k: Number of variables we want our specialised system of multivariate
        polynomials to have.
    :param q: Finite field in which the system of multivariate polynomials are defined.
    :return: Complexity of the Hybrid F5 algorithm
    """
    deg_reg = degree_of_regularity(n - k, m, q)
    if q == 2:
        bi = sum([binomial(n - k, i) for i in range(0, deg_reg + 1)])
    else:
        bi = binomial(n - k - 1 + deg_reg, deg_reg)
    return q**k * bi**2
