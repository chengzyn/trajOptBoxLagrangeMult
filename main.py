from utils import *


def main():
    """Trajectory Optimization of Box Moving Problem Using Lagrange Multiplier
    This code implements the Lagrangian multiplier method for the box moving problem using Lagrange Multipliers,
    where number of collocation points, N, can be adjusted.
    """
    N = 3  # number of collocation points, can adjust this
    M = 3 * N + 2 * N  # number of variables
    assert type(N) is int
    assert N >= 3

    A = np.zeros((M, M))  # matrix containing the coefficients of dL/dz
    final_time = 1
    h = final_time / (N - 1)
    t_vec = [i * h for i in range(N)]

    # A is assembled in 5 parts
    for i in range(0, N):  # dL/dx0 ... dL/dx4
        A[i, 3 * N + i] = -1
        if i <= (N - 3):  # exclude last and last second pt
            A[i, 3 * N + i + 1] = 1
    for i in range(N, 2 * N):  # dL/dv0 ... dL/dv4
        A[i, 3 * N + i] = -1
        if i <= (2 * N - 3):  # exclude last and last second pt
            A[i, 3 * N + i + 1] = 1
        if i == 2 * N - 1:
            pass
        elif i == N:
            A[i, 2 * N + i + 1] = 0.5 * h
        elif i == (2 * N - 2):
            A[i, 2 * N + i] = 0.5 * h
        else:
            A[i, 2 * N + i] = 0.5 * h
            A[i, 2 * N + i + 1] = 0.5 * h
    for i in range(2 * N, 3 * N):  # dL/du0 ... dL/du4
        if i == 2 * N:
            A[i, 2 * i + 1] = 0.5 * h
            A[i, i] = h
        elif i == (3 * N - 1):
            A[i, i] = h
        elif i == (3 * N - 2):
            A[i, 2 * N + i] = 0.5 * h
            A[i, i] = 2 * h
        else:
            A[i, 2 * N + i] = 0.5 * h
            A[i, 2 * N + i + 1] = 0.5 * h
            A[i, i] = 2 * h
    for i in range(3 * N, 4 * N):  # dL/dLambda0 ... dL/dLambda4
        A[i, i - 3 * N] = -1
        if i == 3 * N or i == 4 * N - 1:
            pass
        else:
            A[i, i - 3 * N - 1] = 1
            A[i, i - 2 * N] = 0.5 * h
            A[i, i - 2 * N - 1] = 0.5 * h
    for i in range(4 * N, 5 * N):  # dL/dLambda5 ... dL/dLambda9
        A[i, i - 3 * N] = -1
        if i == 4 * N or i == 5 * N - 1:
            pass
        else:
            A[i, i - 3 * N - 1] = 1
            A[i, i - 2 * N] = 0.5 * h
            A[i, i - 2 * N - 1] = 0.5 * h

    # RHS of Ax=b
    b = np.zeros((M, 1))
    b[M - N - 1, 0] = -1

    # numerical solution
    x = np.linalg.solve(A, b)
    print("computed x:", np.transpose(x[0:N]), "\nanalytical x:", x_analytic(t_vec))
    print("computed u:", np.transpose(x[2 * N:3 * N]), "\nanalytical x:", u_analytic(t_vec))

    # Print out A, add in the indices
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
    main()
