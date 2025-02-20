import numpy as np
from numba import njit


@njit
def legendre_p(n, x):
    """ "
    Returns the value of the Legendre polynomial P_n(x).

    Parameters:
    n (int): The degree of the polynomial.
    x (float): The point at which to evaluate the polynomial. Must satisfy -1 <= x <= 1.
    """
    if n == 0:
        return 1.0
    if n == 1:
        return x

    p_prev2 = 1.0
    p_prev1 = x
    for n in range(2, n + 1):
        p = ((2 * n - 1) * x * p_prev1 - (n - 1) * p_prev2) / n
        p_prev2 = p_prev1
        p_prev1 = p
    return p_prev1


@njit
def assoc_legendre_p(n, m, x):
    """
    Returns the value of the associated Legendre polynomial P_n^m(x).

    Parameters:
    n (int): The degree of the polynomial. Must be greater than or equal to m.
    m (int): The order of the polynomial. Must satisfy 0 <= m <= n.
    x (float): The point at which to evaluate the polynomial. Must satisfy -1 <= x <= 1.
    """
    if m < 0 or m > n:
        return 0.0

    if m == n:
        fact = 1.0
        for k in range(n):
            fact *= -(2 * k + 1)
        return fact * (1 - x * x) ** (n / 2)

    pmm = 1.0
    fact = 1.0
    for k in range(m):
        pmm *= -(2 * k + 1) * (1 - x * x) ** (0.5)

    if n == m:
        return pmm

    pmmp1 = x * (2 * m + 1) * pmm

    if n == m + 1:
        return pmmp1

    for n in range(m + 2, n + 1):
        pmn = (x * (2 * n - 1) * pmmp1 - (n + m - 1) * pmm) / (n - m)
        pmm = pmmp1
        pmmp1 = pmn

    return pmmp1


@njit
def polevl(x: float, coef: np.ndarray, N: int) -> float:
    ans = coef[0]
    i = N
    p = 1  # Indexing in Python instead of pointer

    while i:
        ans = ans * x + coef[p]
        p += 1
        i -= 1

    return ans


@njit
def ellipe(m):
    """
    Compute the complete elliptic integral of the second kind.
    Adapted from Cephes Math Library by Stephen L. Moshier
    MIT License

    Parameters:
        m (float): Input parameter (0 <= x <= 1)

    Returns:
        float: Computed value of E(m)
    """

    P = np.array(
        [
            1.53552577301013293365e-4,
            2.50888492163602060990e-3,
            8.68786816565889628429e-3,
            1.07350949056076193403e-2,
            7.77395492516787092951e-3,
            7.58395289413514708519e-3,
            1.15688436810574127319e-2,
            2.18317996015557253103e-2,
            5.68051945617860553470e-2,
            4.43147180560990850618e-1,
            1.00000000000000000299e0,
        ]
    )

    Q = np.array(
        [
            3.27954898576485872656e-5,
            1.00962792679356715133e-3,
            6.50609489976927491433e-3,
            1.68862163993311317300e-2,
            2.61769742454493659583e-2,
            3.34833904888224918614e-2,
            4.27180926518931511717e-2,
            5.85936634471101055642e-2,
            9.37499997197644278445e-2,
            2.49999999999888314361e-1,
        ]
    )

    x = 1.0 - m
    if x <= 0.0 or x > 1.0:
        if x == 0.0:
            return 1.0
        return 0.0  # Domain error handling

    return polevl(x, P, 10) - np.log(x) * (x * polevl(x, Q, 9))


@njit
def ellipk(m):
    """
    Compute the complete elliptic integral of the second kind.
    Adapted from Cephes Math Library by Stephen L. Moshier
    MIT License

    Parameters:
        m (float): Input parameter (0 <= m <= 1)

    Returns:
        float: Computed value of K(m)
    """

    MACHEP = 1.11022302462515654042e-16

    P = np.array(
        [
            1.37982864606273237150e-4,
            2.28025724005875567385e-3,
            7.97404013220415179367e-3,
            9.85821379021226008714e-3,
            6.87489687449949877925e-3,
            6.18901033637687613229e-3,
            8.79078273952743772254e-3,
            1.49380448916805252718e-2,
            3.08851465246711995998e-2,
            9.65735902811690126535e-2,
            1.38629436111989062502e0,
        ]
    )

    Q = np.array(
        [
            2.94078955048598507511e-5,
            9.14184723865917226571e-4,
            5.94058303753167793257e-3,
            1.54850516649762399335e-2,
            2.39089602715924892727e-2,
            3.01204715227604046988e-2,
            3.73774314173823228969e-2,
            4.88280347570998239232e-2,
            7.03124996963957469739e-2,
            1.24999999999870820058e-1,
            4.99999999999999999821e-1,
        ]
    )

    C1 = 1.3862943611198906188e0  # log(4)

    x = 1.0 - m

    if x < 0.0 or x > 1.0:
        return 0.0

    if x > MACHEP:
        return polevl(x, P, 10) - np.log(x) * polevl(x, Q, 10)
    else:
        if x == 0.0:
            return np.inf
        else:
            return C1 - 0.5 * np.log(x)