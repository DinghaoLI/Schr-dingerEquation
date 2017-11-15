import numpy as np
import matplotlib.pyplot as plt
M = np.loadtxt("plot.txt")
x=M[0]
y=M[1]
plt.plot(x,y,'o')
plt.xlabel('z')
plt.ylabel('Psi')
plt.title("Fonction d'onde pour une equation HO 1D n=0")
plt.show()