from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from sys import exit

r = np.loadtxt("final_config.dat")
r0 = np.loadtxt("initial_config.dat")

ry = r[:,1]
rz = r[:,2]

r0y = r0[:,1]
r0z = r0[:,2]

N = r.shape[0]
for i in range(N):
	plt.plot([r0y[i],ry[i]],[r0z[i],rz[i]])
	plt.scatter(r0y[i],r0z[i],c='red')
	plt.scatter(ry[i],rz[i],c='blue')


plt.show()
exit()
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

for i in range(r.shape[0]):
	ax.scatter(r[i,0],r[i,1],r[i,2])

ax.set_xlim([-.1,10.1])
ax.set_ylim([-.1,10.1])
ax.set_zlim([-.1,10.1])
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

plt.show()



