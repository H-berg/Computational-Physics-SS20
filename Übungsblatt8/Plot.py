import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
plt.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
plt.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
rc('text', usetex=True)


t, x, y, E_kin, E_pot, T = np.genfromtxt("Aequilibrierung.txt", unpack=True)


plt.plot(t, E_kin, label = "kinetische Energie")
plt.plot(t, E_pot, label = "potentielle Energie")

plt.plot(t, E_kin+E_pot, label = 'Gesamtenergie')
plt.legend(loc='best')
plt.ylabel("E")
plt.xlabel("Zeit")
plt.xlim(0, t[-1])
plt.grid()
plt.savefig('Energie.pdf')

plt.clf()

plt.plot(t, T, label = "Temperatur")
plt.legend(loc='best')
plt.ylabel("T")
plt.xlabel("Zeit")
plt.xlim(0, t[-1])
plt.grid()
plt.savefig('Temperatur.pdf')

plt.clf()

plt.plot(t, x, label = "Schwerpunktsgeschwindigkeit x-Richtung")
plt.plot(t, y, label = "Schwerpunktsgeschwindigkeit y-Richtung")
plt.legend(loc='best')
plt.ylabel("v")
plt.xlabel("Zeit")
plt.xlim(0, t[-1])
plt.grid()
plt.savefig('Schwerpunkt.pdf')

plt.clf()

for T_i in [0.01, 1, 100]: #100 erstmal weglassen
    t, x, y, E_kin, E_pot, T = np.genfromtxt("Aequilibrierung"+str(T_i)+".txt", unpack=True)


    plt.plot(t, E_kin, label = "kinetische Energie")
    plt.plot(t, E_pot, label = "potentielle Energie")

    plt.plot(t, E_kin+E_pot, label = 'Gesamtenergie')
    plt.legend(loc='best')
    plt.ylabel("E")
    plt.xlabel("Zeit")
    plt.xlim(0, t[-1])
    plt.grid()
    plt.savefig("Energie"+str(T_i)+".pdf")

    plt.clf()

    plt.plot(t, T, label = "Temperatur")
    plt.legend(loc='best')
    plt.ylabel("T")
    plt.xlabel("Zeit")
    plt.xlim(0, t[-1])
    plt.grid()
    plt.savefig("Temperatur"+str(T_i)+".pdf")

    plt.clf()

    plt.plot(t, x, label = "Schwerpunktsgeschwindigkeit x-Richtung")
    plt.plot(t, y, label = "Schwerpunktsgeschwindigkeit y-Richtung")
    plt.legend(loc='best')
    plt.ylabel("v")
    plt.xlabel("Zeit")
    plt.xlim(0, t[-1])
    plt.grid()
    plt.savefig("Schwerpunkt"+str(T_i)+".pdf")

    plt.clf()

    r, g = np.genfromtxt("Messung"+str(T_i)+".txt", unpack=True)


    plt.plot(r, g, 'bx', label = "Paarkorrelation")
    plt.plot(r, g, 'b')
    plt.legend(loc='best')
    plt.ylabel("g")
    plt.xlabel("r")
    plt.xlim(0, r[-1])
    plt.grid()
    plt.savefig("Korrelation"+str(T_i)+".pdf")

    plt.clf()



