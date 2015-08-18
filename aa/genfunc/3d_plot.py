import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

g = np.genfromtxt("const_energy_2503.dat")
g = g[(g.T[4]<3000)*(g.T[6]<3000)*(g.T[5]<3000)]

lloops = g[g.T[9]>0.5]
sloops = g[g.T[7]>0.5]
boxs = g[(g.T[7]<0.5)*(g.T[8]<0.5)*(g.T[9]<0.5)]

fig = plt.figure(figsize=[6.64,3.32])
ax = fig.add_subplot(1,2, 1, projection='3d')
plt.subplots_adjust(hspace=0.05)
p1 = Line2D([5,6],[7,8],color='k',marker='o',linestyle="none")
p2 = Line2D([5,6],[7,8],color='r',marker='x',linestyle="none")
p3 = Line2D([5,6],[7,8],color='b',markeredgewidth=0.,marker='^',linestyle="none")
leg = plt.legend([p1,p2,p3], ['Short-axis loops','Long-axis loops','Boxes'],loc='upper center',bbox_to_anchor=(0.5,1.1),ncol=3, numpoints = 1)
leg.draw_frame(False)

ax.view_init(azim=47,elev=25)
ax.scatter3D(2.*lloops.T[4],np.abs(lloops.T[5]),np.abs(lloops.T[6]),marker='o', lw = 0,s=5,color='k',label='Short-axis loops',alpha=0.7)
ax.scatter3D(2.*sloops.T[4],np.abs(sloops.T[5]),np.abs(sloops.T[6]),marker='x', lw = 0.5,s=4,color='r',label='Long-axis loops',alpha=0.7)
ax.scatter3D(boxs.T[4],boxs.T[5],boxs.T[6],marker='^', lw = 0,s=6,color='b',label='Boxes',alpha=0.7)
ax.set_xlabel(r"$J_1'/{\rm kpc\,km\,s}^{-1}$")
ax.xaxis.set_major_locator(MaxNLocator(5))
ax.yaxis.set_major_locator(MaxNLocator(5))
ax.zaxis.set_major_locator(MaxNLocator(5))
ax.set_xlim(0,3000)
ax.set_ylabel(r"$J_2'/{\rm kpc\,km\,s}^{-1}$")
ax.set_ylim(0,3000)
ax.set_zlabel(r"$J_3'/{\rm kpc\,km\,s}^{-1}$")
ax.set_zlim(0,3000)
ax.xaxis._axinfo['ticklabel']['space_factor'] =1.1
ax.xaxis._axinfo['label']['space_factor'] = 2.4
ax.yaxis._axinfo['ticklabel']['space_factor'] =1.1
ax.yaxis._axinfo['label']['space_factor'] = 2.4
ax.zaxis._axinfo['ticklabel']['space_factor'] =1.
ax.zaxis._axinfo['label']['space_factor'] = 2.2


ax = fig.add_subplot(1,2, 2, projection='3d')
ax.view_init(azim=156,elev=11)
ax.scatter3D(2.*lloops.T[4],np.abs(lloops.T[5]),np.abs(lloops.T[6]),marker='o', lw = 0,s=5,color='k',alpha=0.7)
ax.scatter3D(2.*sloops.T[4],np.abs(sloops.T[5]),np.abs(sloops.T[6]),marker='x', lw = 0.5,s=4,color='r',alpha=0.7)
ax.scatter3D(boxs.T[4],boxs.T[5],boxs.T[6],marker='^', lw = 0,s=6,color='b',alpha=0.7)
ax.xaxis.set_major_locator(MaxNLocator(5))
ax.yaxis.set_major_locator(MaxNLocator(5))
ax.zaxis.set_major_locator(MaxNLocator(5))
ax.set_xlabel(r"$J_1'/{\rm kpc\,km\,s}^{-1}$")
ax.set_xlim(0,3000)
ax.set_ylabel(r"$J_2'/{\rm kpc\,km\,s}^{-1}$")
ax.set_ylim(0,3000)
ax.set_zlabel(r"$J_3'/{\rm kpc\,km\,s}^{-1}$")
ax.set_zlim(0,3000)
ax.xaxis._axinfo['ticklabel']['space_factor'] = 0.9
ax.xaxis._axinfo['label']['space_factor'] = 2.2
ax.yaxis._axinfo['ticklabel']['space_factor'] =0.8
ax.yaxis._axinfo['label']['space_factor'] = 2.
ax.zaxis._axinfo['ticklabel']['space_factor'] =0.8
ax.zaxis._axinfo['label']['space_factor'] = 1.9
plt.savefig('3d_action_diagram_double.pdf',bbox_inches='tight')
