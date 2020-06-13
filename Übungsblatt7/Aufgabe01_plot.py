import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

plt.rc('text', usetex=True)
x1a,y1a,z1a,t1a, E1 = np.genfromtxt('teil_a1.txt', unpack=True)
x2a,y2a,z2a,t2a, E2 = np.genfromtxt('teil_a2.txt', unpack=True)
E, t = np.genfromtxt('teil_c.txt', unpack=True)

fig1 = plt.figure(1)
fig2 = plt.figure(2)
fig3 = plt.figure(3)
ax1 = fig1.gca(projection='3d')
ax2 = fig2.gca(projection='3d')

p1 = ax1.scatter3D(x1a, y1a, z1a, c=t1a, cmap='viridis')
cbar = fig1.colorbar(p1)
cbar.set_label(r'Zeitschritte', rotation=270, labelpad=10)
ax1.plot(x1a, y1a, z1a) #label=r'$\vec{r}(0) = (1,2,3)^T$ und $\vec{v}(0) = \vec{0}$')
ax1.view_init(70, 30)
plt.xlabel('x')
plt.ylabel('y')
#ax1.legend()
fig1.savefig("Abbildungen/A1_1a_1.pdf")

p2 = ax2.scatter3D(x2a, y2a, z2a, c=t2a, cmap='viridis')
cbar = fig2.colorbar(p2)
cbar.set_label(r'Zeitschritte', rotation=270, labelpad=10)
ax2.plot(x2a, y2a, z2a)# label=r'$\vec{r}(0) = (1,2,3)^T$ und $\vec{v}(0) = (1,1,1)^T$')
ax2.view_init(70, 30)
plt.xlabel('x')
plt.ylabel('y')
#   ax2.legend()
fig2.savefig("Abbildungen/A1_1a_2.pdf")


plt.plot(t1a, E1-E1[0], "bx", label="Energie")
plt.xlabel('Zeit / t')
plt.ylabel('Energie / E')
fig3.savefig("Abbildungen/A1_1c.pdf")
fig3.legend()
plt.clf()
