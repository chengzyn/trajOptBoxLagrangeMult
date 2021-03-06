from utils import *


def manual_3_points():
    """Trajectory Optimization of Box Moving Problem Using Lagrange Multiplier on 3 collocation points
    This code implements the Lagrangian multiplier method for the box moving problem using Lagrange Multipliers,
    where number of collocation points = 3. Coefficients are manually coded in.
    """
    A = np.zeros((15, 15))
    h = 0.5
    t_vec = [i * 0.5 for i in range(3)]

    A[0, 9] = -1
    A[0, 10] = 1

    A[1, 10] = -1

    A[2, 11] = -1

    A[3, 10] = 0.5 * h
    A[3, 12] = -1
    A[3, 13] = 1

    A[4, 10] = 0.5 * h
    A[4, 13] = -1

    A[5, 14] = -1

    A[6, 6] = h
    A[6, 13] = 0.5 * h

    A[7, 7] = 2 * h
    A[7, 13] = 0.5 * h

    A[8, 8] = h

    A[9, 0] = -1

    A[10, 1] = -1
    A[10, 0] = 1
    A[10, 4] = 0.5 * h
    A[10, 3] = 0.5 * h

    A[11, 2] = -1

    A[12, 3] = -1

    A[13, 4] = -1
    A[13, 3] = 1
    A[13, 7] = 0.5 * h
    A[13, 6] = 0.5 * h

    A[14, 5] = -1

    b = np.zeros((15, 1))
    b[11, 0] = -1

    x = np.linalg.solve(A, b)
    print("computed x:", np.transpose(x[0:3]), "\nanalytical x:", x_analytic(t_vec))
    print("computed u:", np.transpose(x[6:9]), "\nanalytical x:", u_analytic(t_vec))

    # Print out A
    N = 3  # number of collocation points
    M = 3 * N + 2 * N  # number of variables
    r_idx = np.arange(0, M)  # this gives a rank 0 array
    c_idx = np.arange(0, M + 1)
    r_idx = np.reshape(r_idx, (-1, 1))  # this gives a rank 1 column array
    c_idx = np.reshape(c_idx, (1, -1))  # this gives a rank 1 row array
    c_idx = c_idx - 1
    c_idx[0, 0] = 0
    A1 = np.concatenate((r_idx, A), axis=1)
    A2 = np.concatenate((c_idx, A1), axis=0)
    mat_print(A2)


if __name__ == '__main__':
    manual_3_points()
