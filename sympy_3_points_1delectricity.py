from sympy import *


def sympy_3_points_added_constraints():
    # declare all the variables symbolically
    V0, V1, V2 = symbols('V0 V1 V2')
    I0, I1, I2 = symbols('I0 I1 I2')
    R0, R1, R2 = symbols('R0 R1 R2')
    l0, l1, l2, l3, l4, l5 = symbols('l0 l1 l2 l3 l4 l5')

    # timestep
    h = 0.5

    # objective function
    f = 0.5 * h * R0 ** 2 + h * R1 ** 2 + 0.5 * h * R2 ** 2

    # constraint functions
    g0 = V0 - 1
    g1 = V1 - V0 - 0.5 * h * (I1 * R1 + I0 * R0)
    g2 = V2 - V1 - 0.5 * h * (I2 * R2 + I1 * R1)
    g3 = V2
    g4 = I0 - I1
    g5 = I1 - I2

    # Lagrangian function
    L = f - l0 * g0 - l1 * g1 - l2 * g2 - l3 * g3 - l4 * g4 - l5 * g5

    # differentials
    dLdV0 = diff(L, V0)
    dLdV1 = diff(L, V1)
    dLdV2 = diff(L, V2)
    dLdI0 = diff(L, I0)
    dLdI1 = diff(L, I1)
    dLdI2 = diff(L, I2)
    dLdR0 = diff(L, R0)
    dLdR1 = diff(L, R1)
    dLdR2 = diff(L, R2)
    dLdl0 = diff(L, l0)
    dLdl1 = diff(L, l1)
    dLdl2 = diff(L, l2)
    dLdl3 = diff(L, l3)
    dLdl4 = diff(L, l4)
    dLdl5 = diff(L, l5)

    # solve
    print(solve(
        [dLdV0, dLdV1, dLdV2, dLdI0, dLdI1, dLdI2, dLdR0, dLdR1, dLdR2, dLdl0, dLdl1, dLdl2, dLdl3, dLdl4,
         dLdl5], [V0, V1, V2, I0, I1, I2, R0, R1, R2, l0, l1, l2, l3, l4, l5]))


if __name__ == '__main__':
    sympy_3_points_added_constraints()
