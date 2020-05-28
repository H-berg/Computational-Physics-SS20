import matplotlib.pyplot as plt
import numpy as np
import math

k, a = np.genfromtxt("fft_1.txt", unpack=True)
x = np.genfromtxt("fft_2.txt", unpack=True)

def f_2(x):
    return (1j)/(np.pi*x)*(((-1)**x)-1)

y = np.linspace(-127,127,128) #ungerade k
z = np.linspace(-126, 128, 128) #gerade k


plt.figure(1, (7, 5))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(k, a, 'g+', label='FFT 1')
plt.xlim(-210, 210)
plt.ylim(-0.1, 0.5)
plt.legend(loc='best', fontsize='small')
plt.xlabel("k")
plt.ylabel("c")
plt.grid()
plt.savefig('a_1_1.jpg')
plt.clf()

plt.figure(1, (7, 5))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(k, x, 'r+', label='FFT 2')
plt.xlim(-210, 210)
plt.ylim(-0.1, 1)
plt.legend(loc='best', fontsize='small')
plt.xlabel("k")
plt.ylabel("c")
plt.grid()
plt.savefig('a_1_2_1.jpg')
plt.clf()

plt.figure(1, (7, 5))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(k, abs(f_2(y)), 'g+', label='Analytische Berechnung')
plt.plot(k, abs(f_2(z)), 'g+')
#plt.plot(k, abs(f_2(m)), 'g+', label='FFT 2')
#plt.plot(k, abs(f_2(n)), 'g+', label='FFT 2')
plt.xlim(-210, 210)
plt.ylim(-0.1, 1)
plt.legend(loc='best', fontsize='small')
plt.xlabel("k")
plt.ylabel("c")
plt.grid()
plt.savefig('a_1_2_2.jpg')
