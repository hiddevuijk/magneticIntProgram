import numpy as np
import matplotlib.pyplot as plt
from sys import exit

bins = np.loadtxt("bins.dat")
rho = np.loadtxt("rho.dat",delimiter=',')
plt.pcolormesh(bins,bins,rho)
plt.colorbar()
plt.show()

exit()

py = np.loadtxt("py.dat",delimiter=',')
px = np.loadtxt("px.dat",delimiter=',')


plt.plot(bins,py[:,0],color="blue")
plt.plot(bins,py[:,1],color="red")
plt.plot(bins,py[:,2],color="green")


plt.show()





