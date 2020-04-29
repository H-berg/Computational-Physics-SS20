import numpy as np
import matplotlib.pyplot as plt

m, n = np.genfromtxt('data2c.txt', unpack=True)
x, y = np.genfromtxt('data2a.txt', unpack=True)

x_0 = np.linspace(-7, 13, 100)
plt.plot(x_0, x_0*m+n, label = 'Ausgleichsgerade')
plt.plot(x, y, '+r' , label = 'Datenpunkte')
plt.ylabel('y')
plt.xlabel('x')
plt.legend(loc="best")
plt.savefig('plot2d.pdf')
