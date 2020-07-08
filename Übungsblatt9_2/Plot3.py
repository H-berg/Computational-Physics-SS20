import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
plt.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
plt.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
rc('text', usetex=True)

x, y = np.genfromtxt("Initialisierung.txt", unpack=True)

plt.plot(x, y, color = "black", linestyle ='-', linewidth = 0.1)
plt.plot(x, y, 'r.', label = "Ortsvektoren")
plt.plot(0, 0, 'bo', label = "Start")
plt.plot(0.2, 0, 'go', label = "Ende")

plt.legend(loc='best')
plt.ylabel("y")
plt.xlabel("x")
plt.savefig('Initialisierung.pdf')

plt.clf()

x, y = np.genfromtxt("Startkonfig.txt", unpack=True)

plt.plot(x, y, color = "black", linestyle ='-', linewidth = 0.1)
plt.plot(x, y, 'r.', label = "Ortsvektoren")
plt.plot(x[0], y[0], 'bo', label = "Start")
plt.plot(x[89], y[89], 'go', label = "Ende")

plt.legend(loc='best')
plt.ylabel("y")
plt.xlabel("x")
plt.savefig('Startkonfig.pdf')

plt.clf()


for d in [0.9, 0.99, 0.999]:
    for S in [10, 100, 1000, 10000]:
        x, y = np.genfromtxt("Strecke_d_"+str(d)+"_S_"+str(S)+".txt", unpack=True)

        plt.plot(x, y, color = "black", linestyle ='-', linewidth = 0.1)
        plt.plot(x, y, 'r.', label = "Ortsvektoren")
        plt.plot(x[0], y[0], 'bo', label = "Start")
        plt.plot(x[89], y[89], 'go', label = "Ende")

        plt.legend(loc='best')
        plt.ylabel("y")
        plt.xlabel("x")
        plt.savefig("Strecke_d_"+str(d)+"_S_"+str(S)+".pdf")

        plt.clf()
