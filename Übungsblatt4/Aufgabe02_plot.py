import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.rc('text', usetex=True)

a_out = pd.read_csv('out.txt', skiprows=1, decimal='.', delimiter=';')
a_in = pd.read_csv('in.txt', skiprows=1, decimal='.', delimiter=';')

x = np.linspace(1, 8, 1000)
plt.plot(a_out.x, a_out.phi_a, 'r.', label=r"$\phi(x>a)$")
plt.plot(x, 8/x, label=r'$\frac{1}{x}$')
plt.legend(loc='best')
plt.xlabel(r'$\frac{x}{a}$')
plt.ylabel(r'$\Phi(\frac{x}{a})$')
plt.savefig("Abbildungen/out_a.pdf")
plt.clf()

plt.plot(a_in.x, a_in.phi_a, 'r.', label=r"$\phi(x\leq a)$")
plt.legend(loc='best')
plt.xlabel(r'$\frac{x}{a}$')
plt.ylabel(r'$\Phi(\frac{x}{a})$')
plt.savefig("Abbildungen/in_a.pdf")
plt.clf()

plt.plot(a_out.x, a_out.phi_b, 'r.', label=r"$\phi(x>a)$")
plt.plot(x, 8/(3*x**2), label=r'$\frac{8}{3x}$')
plt.legend(loc='best')
plt.xlabel(r'$\frac{x}{a}$')
plt.ylabel(r'$\Phi(\frac{x}{a})$')
plt.savefig("Abbildungen/out_b.pdf")
plt.clf()

plt.plot(a_in.x, a_in.phi_b, 'r.', label=r"$\phi(x\leq a)$")
plt.legend(loc='best')
plt.xlabel(r'$\frac{x}{a}$')
plt.ylabel(r'$\Phi(\frac{x}{a})$')
plt.savefig("Abbildungen/in_b.pdf")
plt.clf()
