import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
plt.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
plt.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
rc('text', usetex=True)

x1, x2 , err1 = np.genfromtxt("data_1.txt", unpack=True)

xlist = np.linspace(-1.1, 1.1, 1000)
ylist = np.linspace(-0.1, 1.1, 1000)
X, Y = np.meshgrid(xlist, ylist)
Z = (1-X)**2+100*(Y-X**2)**2
fig,ax=plt.subplots(1,1)
cp = ax.contourf(X, Y, Z)
fig.colorbar(cp) # Add a colorbar to a plot
plt.plot(x1, x2, label = 'Gradientenverfahren')
plt.legend(loc='best')
plt.savefig('plotgrad_a.pdf')

x1, x2, err2 = np.genfromtxt("data_2.txt", unpack=True)

xlist = np.linspace(-5, 3, 1000)
ylist = np.linspace(-2, 5, 1000)
X, Y = np.meshgrid(xlist, ylist)
Z = (1-X)**2+100*(Y-X**2)**2
fig,ax=plt.subplots(1,1)
cp = ax.contourf(X, Y, Z,)
fig.colorbar(cp) # Add a colorbar to a plot
plt.plot(x1, x2, label = 'konjugiertes Gradientenverfahren')
plt.legend(loc='best')
plt.savefig('plotgrad_b.pdf')

plt.clf()

#Plotten des Fehlers
size = np.size(err1)
x_1 = np.linspace(0, size, size)
size = np.size(err2)
x_2 = np.linspace(0, size, size)

plt.plot(x_1, err1, label = "Gradientenverfahren")
plt.plot(x_2, err2, label = "konjugiertes Gradientenverfahren")
plt.legend(loc='best')
plt.ylabel("Fehler")
plt.xlabel("Anzahl Iterationsschritte")
plt.grid()
plt.savefig('err.pdf')
