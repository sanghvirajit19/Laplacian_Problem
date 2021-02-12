import numpy as np
import matplotlib.pyplot as plt
import seaborn as sn
import timeit
import sys
from numpy import linalg as LA

np.set_printoptions(threshold=sys.maxsize)

def Grid_with_BC(n):

    #grid
    f = np.zeros((n, n))

    # Boundary conditions

    # Top
    f[0, :(n // 4)] = np.arange(13, 5, -(13 - 5) / (n // 4))
    f[1, :(n // 4)] = np.arange(13, 5, -(13 - 5) / (n // 4))

    f[:2, (n // 4):(3 * n // 4)] = 5

    f[0, (3 * n // 4):] = np.arange(5, 13, (13 - 5) / (n // 4))
    f[1, (3 * n // 4):] = np.arange(5, 13, (13 - 5) / (n // 4))

    # Bottom
    f[n-2:, :] = 21

    # Left
    f[:(3 * n // 8), 0] = np.arange(13, 40, ((40 - 13) / (3 * n // 8)))
    f[:(3 * n // 8), 1] = np.arange(13, 40, ((40 - 13) / (3 * n // 8)))

    f[(n // 2):, 0] = np.arange(40, 21, -((40 - 21) / (n // 2)))
    f[(n // 2):, 1] = np.arange(40, 21, -((40 - 21) / (n // 2)))

    # Right

    f[:(n // 2), -1] = np.arange(13, 40, ((40 - 13) / (n // 2)))
    f[:(n // 2), -2] = np.arange(13, 40, ((40 - 13) / (n // 2)))

    f[(5 * n // 8):, -1] = np.arange(40, 21, -((40 - 21) / (3 * n // 8)))
    f[(5 * n // 8):, -2] = np.arange(40, 21, -((40 - 21) / (3 * n // 8)))

    # Heater
    f[(3 * n // 8):(n // 2) + 1, :(n // 8 + 1)] = 40

    f[(n // 2):(5 * n // 8) + 1, -(n // 8 + 1):] = 40

    return f

def jacobi_step(T):

    m, n = T.shape

    _T = np.copy(T)

    # iterate over interior

    for i in range(2, m-1):
        for j in range(2, n-1):

            _T[i, j] = (T[i+1, j] + T[i-1, j] + T[i, j-1] + T[i, j+1]) / 4

    return _T

def restriction(T):

    m = T.shape[0] // 2
    n = T.shape[1] // 2

    _T = np.zeros((m, n))
    s = 2                   #Taking only alternate rows and columns

    # Heater should be there at its own place all the time
    _T[(3 * n // 8):(n // 2) + 1, :(n // 8 + 1)] = 40
    _T[(n // 2):(5 * n // 8) + 1, -(n // 8 + 1):] = 40

    for i in range(m):
        for j in range(n):
            _T[i, j] = T[i * s, j * s]

    return _T

def prolongation(T):

    m = T.shape[0] * 2
    n = T.shape[1] * 2

    _T = Grid_with_BC(m)

    # Heater should be there at its own place all the time
    _T[(3 * n // 8):(n // 2) + 1, :(n // 8 + 1)] = 40
    _T[(n // 2):(5 * n // 8) + 1, -(n // 8 + 1):] = 40

    for i in range(2, m-2):
        for j in range(2, n-2):

            _T[i, j] = (T[int(np.floor((i+1)/2)), int(np.floor((j+1)/2))]
                       + T[int(np.ceil((i+1)/2)), int(np.floor((j+1)/2))]
                       + T[int(np.floor((i+1)/2)), int(np.ceil((j+1)/2))]
                       + T[int(np.ceil((i+1)/2)), int(np.ceil((j+1)/2))]) / 4

    return _T

def simulation(T, epsilon=None, num_steps=None, residual_plot=False, runtime=False):

    start = timeit.default_timer()
    global _T, residual, iteration, residual_list
    m, n = T.shape

    iteration_list = []
    l2_list = []
    residual = 1

    for i in range(num_steps):

        i += 1
        # Heater should be there at its own place all the time
        T[(3 * n // 8):(n // 2) + 1, :(n // 8 + 1)] = 40
        T[(n // 2):(5 * n // 8) + 1, -(n // 8 + 1):] = 40

        _T = jacobi_step(T)
        error = _T - T
        T = _T

        l2_norm = LA.norm(error)
        l2_list.append(l2_norm)

        iteration_list.append(i)

        if residual < np.log10(epsilon):
            print("Convergence Criteria satisfied")
            print("Solution converged in {} iterations".format(i))
            break

        l2_array = np.asarray(l2_list)
        maximum_value = l2_array.max()
        resi_drop = l2_array / maximum_value

        residual_list = np.log(resi_drop).reshape(-1, 1)
        iteration = np.asarray(iteration_list).reshape(-1, 1)

        residual = residual_list[-1]

        if i == num_steps-1:
            print("Iteration criteria satisfied")
            print("Solution converged in {} iterations".format(i))

    if residual_plot == True:
        plt.plot(iteration, residual_list)
        plt.xlabel("Iterations")
        plt.ylabel("Residual drop")
        plt.show()

    end = timeit.default_timer()

    if runtime == True:
        print("runtime: {} s".format(float(round(end - start, 3))))

    return _T

def single_level(grid, runtime=False):

    start = timeit.default_timer()
    # Making the grid with BC
    A = Grid_with_BC(grid)

    # Visualizing the problem
    sn.heatmap(A)
    plt.show()

    # Running the simulation
    B = simulation(A, epsilon=1e-5, num_steps=10000, residual_plot=True)

    # Approximate solution
    sn.heatmap(B)
    plt.show()

    end = timeit.default_timer()

    if runtime == True:
        print("runtime: {} s".format(float(round(end - start, 3))))

    return B

def multi_level(grid, runtime=False):

    start_ = timeit.default_timer()
    print("Simulation Running...")
    A = Grid_with_BC(grid)

    # Approximate solution
    sn.heatmap(A)
    plt.show()

    B = simulation(A, epsilon=1e-5, num_steps=5000, residual_plot=True, runtime=False)

    # Approximate solution
    sn.heatmap(B)
    plt.show()

    # Downlevel
    C = restriction(B)
    print("downleveled")

    # Approximate solution
    sn.heatmap(C)
    plt.show()

    D = simulation(C, epsilon=1e-5, num_steps=5000, residual_plot=True, runtime=False)

    # Approximate solution
    sn.heatmap(D)
    plt.show()

    # Uplevel
    E = prolongation(D)
    print("upleveled")

    # Approximate solution
    sn.heatmap(E)
    plt.show()

    F = simulation(E, epsilon=1e-5, num_steps=5000, residual_plot=True, runtime=False)

    # Approximate solution
    sn.heatmap(F)
    plt.show()

    end_ = timeit.default_timer()

    if runtime == True:
        print("runtime: {} s".format(float(round(end_ - start_, 3))))

    return E

solution_A = single_level(256, runtime=True)

solution_B = multi_level(256, runtime=True)
