import numpy as np
import matplotlib.pyplot as plt

x = np.loadtxt("x.dat")
y = np.loadtxt("y.dat")

plt.plot(x,y)
plt.show()


