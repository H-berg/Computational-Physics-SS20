import matplotlib.pyplot as plt
import numpy as np

a_direkt, a_verb, rel_err_a = np.genfromtxt("a_2_a.txt", unpack=True)

plt.figure(1, (7, 5))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(a_direkt, 'r-', label='Direkte Rechnung')
plt.plot(a_verb, 'g-', label='Verbesserte Rechnung')
plt.legend(loc='best', fontsize='x-small')
plt.savefig('a_2_a.pdf')
plt.clf()

b_direkt, b_verb, rel_err_b = np.genfromtxt("a_2_b.txt", unpack=True)

plt.figure(1, (7, 5))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(b_direkt, 'r-', label='Direkte Rechnung')
plt.plot(b_verb, 'g-', label='Verbesserte Rechnung')
plt.legend(loc='best', fontsize='x-small')
plt.savefig('a_2_b.pdf')
plt.clf()

c_direkt, c_verb, rel_err_c = np.genfromtxt("a_2_c.txt", unpack=True)

plt.figure(1, (7, 5))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(c_direkt, 'r-', label='Direkte Rechnung')
plt.plot(c_verb, 'g-', label='Verbesserte Rechnung')
plt.legend(loc='best', fontsize='x-small')
plt.savefig('a_2_c.pdf')
plt.clf()

plt.figure(1, (7, 5))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.subplot(3,1,1)
plt.plot(rel_err_a, 'r-', label='Aufgabenteil a')
plt.legend(loc='best', fontsize='x-small')
plt.subplot(3,1,2)
plt.plot(rel_err_b, 'g-', label='Aufgabenteil b')
plt.legend(loc='best', fontsize='x-small')
plt.subplot(3,1,3)
plt.plot(rel_err_c, 'b-', label='Aufgabenteil c')
plt.legend(loc='best', fontsize='x-small')
plt.savefig('a_2_rel_err.pdf')
