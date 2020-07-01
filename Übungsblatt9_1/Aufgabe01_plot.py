import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

plt.rc('text', usetex=True)
H, m = np.genfromtxt('A1.txt', unpack=True)

m_max = max(abs(m))

# normieren
m /= m_max

plt.plot(H, m,  linewidth=6, color='lightcoral', label="Numerisch")
plt.plot(H, np.tanh(H) ,color='#4b0082', label="Analytisch")
plt.axis([-5.3, 5.3, -1.1, 1.1])
plt.xticks(np.linspace(-5,5,6))
plt.xlabel("Magnetfeldstärke H")
plt.ylabel("Magnetisierung m")
plt.legend()
#plt.show()
plt.savefig("Abbildungen/tanh.pdf")
plt.clf()
