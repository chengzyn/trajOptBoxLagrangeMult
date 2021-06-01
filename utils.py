import numpy as np


def u_analytic(t):
    """Computes the analytic solution for u (force) at time t."""
    return 6 - 12 * np.array(t)


def x_analytic(t):
    """Computes the analytic solution for x (position) at time t."""
    t = np.array(t)
    return 3 * np.power(t, 2) - 2 * np.power(t, 3)


def mat_print(mat, fmt="g"):
    """Prints a matrix out with equal spacing for visual checking."""
    col_maxes = [max([len(("{:" + fmt + "}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:" + str(col_maxes[i]) + fmt + "}").format(y), end="  ")
        print("")
