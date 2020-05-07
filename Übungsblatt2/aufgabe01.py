import numpy as np
import matplotlib.pyplot as plt

for x in np.array([0,10,20,50]):
    name = "bild_rang_" + str(x)
    A = np.genfromtxt(name + ".txt")
    plt.imshow(A, cmap='gray')
    plt.savefig(name + ".pdf")
