from sympy import *


def sympy_3_points_added_constraints():
    # declare all the variables symbolically
    u0, u1, u2 = symbols('u0 u1 u2')
    x0, x1, x2 = symbols('x0 x1 x2')
    v0, v1, v2 = symbols('v0 v1 v2')
    l0, l1, l2, l3, l4, l5, l6, l7 = symbols('l0 l1 l2 l3 l4 l5 l6 l7')

    # timestep
    h = 0.5

    # objective function
    f = 0.5 * h * u0 ** 2 + h * u1 ** 2 + 0.5 * h * u2 ** 2

    # constraint functions
    g0 = x0
    g1 = x1 - x0 - 0.5 * h * (v1 + v0)
    g2 = x2 - 1
    g3 = v0
    g4 = v1 - v0 - 0.5 * h * (u1 + u0)
    g5 = v2
    g6 = x2 - x1 - 0.5 * h * (v2 + v1)
    g7 = v2 - v1 - 0.5 * h * (u2 + u1)

    # Lagrangian function
    L = f - l0 * g0 - l1 * g1 - l2 * g2 - l3 * g3 - l4 * g4 - l5 * g5 - l6 * g6 - l7 * g7

    # differentials
    dLdx0 = diff(L, x0)
    dLdx1 = diff(L, x1)
    dLdx2 = diff(L, x2)
    dLdv0 = diff(L, v0)
    dLdv1 = diff(L, v1)
    dLdv2 = diff(L, v2)
    dLdu0 = diff(L, u0)
    dLdu1 = diff(L, u1)
    dLdu2 = diff(L, u2)
    dLdl0 = diff(L, l0)
    dLdl1 = diff(L, l1)
    dLdl2 = diff(L, l2)
    dLdl3 = diff(L, l3)
    dLdl4 = diff(L, l4)
    dLdl5 = diff(L, l5)
    dLdl6 = diff(L, l6)
    dLdl7 = diff(L, l7)

    # solve
    print(linsolve(
        [dLdx0, dLdx1, dLdx2, dLdv0, dLdv1, dLdv2, dLdu0, dLdu1, dLdu2, dLdl0, dLdl1, dLdl2, dLdl3, dLdl4, dLdl5, dLdl6,
         dLdl7], (x0, x1, x2, v0, v1, v2, u0, u1, u2, l0, l1, l2, l3, l4, l5, l6, l7)))


if __name__ == '__main__':
    sympy_3_points_added_constraints()
