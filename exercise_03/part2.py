import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

# load data
n, theta, t_tree, m_used_nodes, t_exact, err = np.loadtxt("results.txt", skiprows=1, unpack=True)



# loop over angles
for th in theta:
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    fig.suptitle("$\\theta = {:1.2f}$".format(th))

    ax1.plot(n[theta == th], t_exact[theta == th], label="Exact summation")
    ax1.plot(n[theta == th], t_tree[theta == th], label="Tree algorhitm")
    ax1.set_xscale("log")
    ax1.set_title("Duration")
    ax1.set_ylabel("Time [sec]")
    ax1.legend(loc=2)

    ax2.plot(n[theta == th], err[theta == th])
    ax2.set_title("Mean Error per particle")
    ax2.set_xlabel("Number of particles")
    ax2.legend(loc=4)

    plt.savefig("eval_th_{:d}.png".format(int(th*10)))


# calculation of 10^10 particles
coefficient_ex = t_exact/(n**2)
coefficient_tree = t_tree/(np.log(n)*n)
print(coefficient_tree)
print("1e10 particles")
print("\texact: {:3.10f} * N^2 = {:25f}sec".format(np.mean(coefficient_ex), np.mean(coefficient_ex)*(1e20)))
print("\ttree: {:3.10f} * ln(N)*N = {:25f}sec".format(np.mean(coefficient_tree), np.mean(coefficient_tree)*np.log(1e10)*1e10))
