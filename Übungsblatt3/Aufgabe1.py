import numpy as np
import matplotlib.pyplot as plt

A = np.genfromtxt('Bild.txt')
plt.imshow(A, cmap='gray')
plt.savefig('Bild1b.pdf')

A10 = np.genfromtxt('data1b10.txt')
plt.imshow(A10, cmap='gray')
plt.savefig('Bild1b10.pdf')

A20 = np.genfromtxt('data1b20.txt')
plt.imshow(A20, cmap='gray')
plt.savefig('Bild1b20.pdf')

A50 = np.genfromtxt('data1b50.txt')
plt.imshow(A50, cmap='gray')
plt.savefig('Bild1b50.pdf')
