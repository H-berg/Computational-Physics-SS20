import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

plt.rc('text', usetex=True)
x1a,y1a,z1a, t, E = np.genfromtxt('A2_1a.txt', unpack=True)
x2a,y2a,z2a = np.genfromtxt('A2_2a.txt', unpack=True)
x3a,y3a,z3a = np.genfromtxt('A2_3a.txt', unpack=True)

fig1 = plt.figure(1)
fig2 = plt.figure(2)
fig3 = plt.figure(3)
fig4 = plt.figure(4)
ax1 = fig1.gca(projection='3d')
ax2 = fig2.gca(projection='3d')
ax3 = fig3.gca(projection='3d')

p1 = ax1.scatter3D(x1a, y1a, z1a, c=t, cmap='viridis')
cbar = fig1.colorbar(p1)
cbar.set_label(r'Zeitschritte', rotation=270, labelpad=10)
ax1.plot(x1a, y1a, z1a,)
ax1.view_init(55, 80)
plt.xlabel('x')
plt.ylabel('y')
fig1.savefig("Abbildungen/A2_a1.pdf")

p2 = ax2.scatter3D(x2a, y2a, z2a, c=t, cmap='viridis')
cbar = fig2.colorbar(p2)
cbar.set_label(r'Zeitschritte', rotation=270, labelpad=10)
ax2.plot(x2a, y2a, z2a)
ax2.view_init(55, 80)
plt.xlabel('x')
plt.ylabel('y')
fig2.savefig("Abbildungen/A2_a2.pdf")

p3 = ax3.scatter3D(x3a, y3a, z3a, c=t, cmap='viridis')
cbar = fig3.colorbar(p3)
cbar.set_label(r'Zeitschritte', rotation=270, labelpad=10)
ax3.plot(x3a, y3a, z3a)
ax3.view_init(55, 80)
plt.xlabel('x')
plt.ylabel('y')
fig3.savefig("Abbildungen/A2_a3.pdf")

plt.plot(t, E, label="Energie")
plt.xlabel('Zeit / t')
plt.ylabel('Energie / E')
fig4.savefig("Abbildungen/A2_b.pdf")
plt.clf()
