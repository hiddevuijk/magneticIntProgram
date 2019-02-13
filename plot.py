import numpy as np
import matplotlib.pyplot as plt
from sys import exit

plt.subplot(1,2,1)
presX = np.loadtxt("pressureX.dat")
plt.imshow(presX,interpolation="none")
plt.colorbar()

plt.subplot(1,2,2)
presY = np.loadtxt("pressureY.dat")
plt.imshow(presY,interpolation="none")
plt.colorbar()

plt.show()


rho = np.loadtxt("rho.dat")
bins = np.loadtxt("rho_bins.dat")

fx = np.loadtxt("fx.dat")
fy = np.loadtxt("fy.dat")
fz = np.loadtxt("fz.dat")
print(sum(sum(fy)))
px = np.loadtxt("px.dat")
py = np.loadtxt("py.dat")
pz = np.loadtxt("pz.dat")

plt.subplot(2,3,1)
plt.imshow(rho,interpolation='none')
plt.title("density")
plt.colorbar()

plt.subplot(2,3,2)
plt.imshow(fx,interpolation='none')
plt.title("flux x")

plt.subplot(2,3,3)
plt.imshow(fy,interpolation='none')
plt.title("flux Y")

plt.subplot(2,3,5)
plt.imshow(px,interpolation='none')
plt.title("p x")

plt.subplot(2,3,6)
plt.imshow(py,interpolation='none')
plt.title("p y")

plt.show()


