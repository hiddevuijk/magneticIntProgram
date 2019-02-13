import numpy as np
import matplotlib.pyplot as plt
from sys import exit

rho = np.loadtxt("rho.dat")
bins = np.loadtxt("p_bins.dat")

fx = np.loadtxt("fx.dat")
fy = np.loadtxt("fy.dat")
fz = np.loadtxt("fz.dat")

px = np.loadtxt("px.dat")
py = np.loadtxt("py.dat")
pz = np.loadtxt("pz.dat")

plt.subplot(1,3,1)
rhoY = np.mean(rho,axis=1)
plt.plot(bins,rhoY)

plt.subplot(1,3,2)
fxY = np.mean(fx,axis=1)
plt.plot(bins,fxY,label="fx")

fyY = np.mean(fy,axis=1)
plt.plot(bins,fyY,label="fy")
plt.legend()

plt.subplot(1,3,3)
pxY = np.mean(px,axis=1)
plt.plot(bins,pxY,label="px")

pyY = np.mean(py,axis=1)
plt.plot(bins,pyY,label="py")
plt.legend()

plt.show()


