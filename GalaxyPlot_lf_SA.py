import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


x=[]
y=[]
z=[]

f1=open('xy.dat','r')
for line in f1.readlines():
	line=line.split()
        x.append(float(line[0]))
        y.append(float(line[1]))
        z.append(float(line[2]))
f1.close()

f2=open('galaxy_lf.dat','r')
for line in f2.readlines():
	line=line.split()
	x.append(float(line[1]))
	y.append(float(line[2]))
	z.append(float(line[3]))	
f2.close()

for n in range(0,1):
	name="galaxylf_%d.jpg"%n
	min=(1000*n)
	med=(1000*(n+1)-500)
	max=(1000*(n+1))
        fig=plt.figure(n)
        ax=Axes3D(fig)
        ax.set_xlim(-60,60)
        ax.set_ylim(-60,60)
        ax.set_zlim(-60,60)
        #print x[499]
        ax.scatter(x[min:med],y[min:med],z[min:med],c='r')
        ax.scatter(x[med:max],y[med:max],z[med:max],c='b')
        ax.set_xlabel('X position')
        ax.set_ylabel('Y position')
        ax.set_zlabel('Z position')

        plt.savefig(name)
        plt.close(name)
        plt.show()

