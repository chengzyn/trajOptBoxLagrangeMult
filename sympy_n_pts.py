from sympy import *
from utils import *
import matplotlib.pyplot as plt


def sympy_n_pts(N):
    """Trajectory Optimization of Box Moving Problem Using Lagrange Multiplier
    This code implements the Lagrangian multiplier method for the box moving problem using Lagrange Multipliers,
    where number of collocation points, N, can be adjusted.
    """
    # set up collocation
    assert type(N) is int
    assert N >= 3
    final_time = 1
    h = final_time / (N - 1)
    t_vec = [i * h for i in range(N)]

    # declare symbols
    x = [symbols('x%d' % i) for i in range(N)]
    v = [symbols('v%d' % i) for i in range(N)]
    u = [symbols('u%d' % i) for i in range(N)]
    l = [symbols('l%d' % i) for i in range((N + 1) * 2)]

    # objective function
    f1 = [h * pow(i, 2) for i in u]
    f1[0] *= 0.5
    f1[-1] *= 0.5
    f = sum(f1)

    # constraint functions
    g = list(range(len(l)))
    for i in range(len(l)):
        if i == 0:
            g[i] = x[i]
        elif 0 < i < N:
            g[i] = x[i] - x[i - 1] - 0.5 * h * (v[i] + v[i - 1])
        elif i == N:
            g[i] = x[i - 1] - 1
        elif i == N + 1:
            g[i] = v[0]
        elif (N + 1) < i < (2 * N + 1):
            j = i - (N + 1)
            g[i] = v[j] - v[j - 1] - 0.5 * h * (u[j] + u[j - 1])
        else:
            g[i] = v[N - 1]

    # Lagrangian function
    lg = [i * j for i, j in zip(l, g)]
    L = f - sum(lg)

    # differentials
    z = x + v + u + l
    dL = [diff(L, i) for i in z]
    z_comp = linsolve(dL, tuple(z))
    x_comp = z_comp.args[0][0:N]
    v_comp = z_comp.args[0][N:2 * N]
    u_comp = z_comp.args[0][2 * N:3 * N]

    return x_comp, v_comp, u_comp, x_analytic(t_vec), u_analytic(t_vec), t_vec, N


def nice_plot(xc, vc, uc, xa, ua, t, n):
    """Plots the results"""
    mse_x = mean_squared_error(xa, xc)
    mse_u = mean_squared_error(ua, uc)
    fig, axs = plt.subplots(2)
    axs[0].plot(t, xa, label='x_analytical')
    axs[0].plot(t, xc, 'ro', label=f'x_computed\n MSE = {mse_x}', markersize=2)
    axs[0].legend()
    axs[0].set_title(f'N = {n}')
    axs[0].set(xlabel='t (s)', ylabel='x')
    axs[1].plot(t, ua, 'g', label='u_analytical')
    axs[1].plot(t, uc, 'mo', label=f'u_computed\n MSE = {mse_u}', markersize=2)
    axs[1].legend()
    axs[1].set(xlabel='t (s)', ylabel='u')

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()
    plt.savefig(f'n_is_{n}.png', dpi=300)  # must be before show()
    plt.show()


if __name__ == '__main__':
    nice_plot(*sympy_n_pts(51))
