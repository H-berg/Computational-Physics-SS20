import numpy as np
import matplotlib.pyplot as plt

x, euler, symeuler, analytisch = np.genfromtxt('data3a.txt', unpack=True)

plt.plot(x, euler, label='Euler Verfahren')
plt.plot(x, symeuler, label = 'symmetrisches Euler Verfahren')
plt.plot(x, analytisch, label = 'exp(-x)')
plt.ylabel('y')
plt.xlabel('x')
plt.legend(loc="best")
plt.savefig('plot3a.pdf')
