import numpy as np
import matplotlib.pyplot as plt

xb, eulerb, symeulerb, analytischb = np.genfromtxt('eulerb.txt', unpack=True)

plt.plot(xb, eulerb, label='Euler Verfahren')
plt.plot(xb, symeulerb, label = 'symmetrisches Euler Verfahren')
plt.plot(xb, analytischb, label = 'exp(-x)')
plt.ylabel('y')
plt.xlabel('x')
plt.legend(loc="best")
plt.savefig('eulerb.pdf')
