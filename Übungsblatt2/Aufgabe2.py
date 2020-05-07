import numpy as np
import matplotlib.pyplot as plt

t1, t2, t3 = np.genfromtxt('data2b.txt', unpack=True)
N = ([1, 2, 2**2, 2**3, 2**4, 2**5, 2**6, 2**7, 2**8, 2**9, 2**(10), 2**(11), 2**(12), 2**(13)])

plt.figure(1, (7, 5))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.loglog(N, t1, label = 'Erstellen der Zufallsmatrix')
plt.loglog(N, t2, label = 'LU-Zerlegung')
plt.loglog(N, t3, label = 'Loesen des Gleichungssystems')
plt.loglog(N, t3+t2+t1, label = 'Gesamtlaufzeit')
plt.ylabel('Zeit/s')
plt.xlabel('N')
plt.legend(loc="best")
plt.savefig('plot2b.pdf')

#In einem nicht logarithmischen Plot sieht man nochmal dass die LU Zerlegung die meisten Laufzeit beansprucht
plt.clf()
plt.figure(1, (7, 5))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
x = np.linspace(0, 8000, 5000)
plt.plot(N, t1, label = 'Erstellen der Zufallskomponenten')
plt.plot(N, t2, label = 'LU-Zerlegung')
plt.plot(N, t3, label = 'Loesen des Gleichungssystems')
plt.plot(N, t3+t2+t1, label = 'Gesamtlaufzeit')
plt.ylabel('Zeit/s')
plt.xlabel('N')
plt.legend(loc="best")
plt.savefig('plot2c.pdf')
