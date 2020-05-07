import matplotlib.pyplot as plt
import numpy as np

n, inv, fullpiv, partialpiv= np.genfromtxt("times.txt", unpack=True)

plt.figure(1, (7, 5))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(n, inv, 'b-', label='Berechnung mit inverser Matrix')
plt.plot(n, fullpiv, 'c-', label='Volle Pivotisierung')
plt.plot(n, partialpiv, 'm-', label='Teilweise Pivotisierung')
plt.yscale('log')
plt.xscale('log')
plt.ylabel('t/s')
plt.xlabel('N')
plt.legend(loc='best', fontsize='small')
plt.savefig('a_3.pdf')
