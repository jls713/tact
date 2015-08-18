import numpy as np
import matplotlib.pyplot as plt
from matplotlib import lines
import sys
import seaborn as sns

results = np.genfromtxt(sys.argv[1])
labels = ['Fudge','ItTorus','Genfunc','AvGenfunc','uvOrb','PAA','SAA','FIT']
f,a=plt.subplots(3,1,figsize=(3,5))
plt.subplots_adjust(hspace=0.)
JR = np.mean(results.T[4])
Jz = np.mean(results.T[5])
for i in range(3):
	a[i].set_xlim(0.,1000.)
for i in range(8):
	n=i
	if(i>5):
		n=i-6
	l1, = a[0].plot(results.T[i*2],label=labels[i],lw=.5)
	l2, = a[1].plot(results.T[i*2+1],lw=0.5)
	if(i>4):
		l1.set_dashes((3,2))
		l2.set_dashes((3,2))
	JRspread = np.std(results.T[i*2])
	x,y=np.array([[1050+10*i,1050+10*i],[JR-JRspread,JR+JRspread]])
	line = lines.Line2D(x, y, lw=1.,color=sns.color_palette()[n])
	if(i>4):
		line.set_dashes((3,2))
	line.set_clip_on(False)
	a[0].add_line(line)
	Jzspread = np.std(results.T[i*2+1])
	x,y=np.array([[1050+10*i,1050+10*i],[Jz-Jzspread,Jz+Jzspread]])
	line = lines.Line2D(x, y, lw=1.,color=sns.color_palette()[n])
	if(i>4):
		line.set_dashes((3,2))
	line.set_clip_on(False)
	a[1].add_line(line)

plt.setp(a[0].get_xticklabels(),visible=False)
plt.setp(a[1].get_xticklabels(),visible=False)
R = np.sqrt(results.T[16]*results.T[16]+results.T[17]*results.T[17])
a[2].plot(R,lw=.5,color='k')
a[0].set_xlabel(r'$t$')
a[1].set_xlabel(r'$t$')
a[0].set_ylabel(r'$J_R$')
a[1].set_ylabel(r'$J_z$')
a[2].set_ylabel(r'$R$')
a[2].set_xlabel(r'$t$')

a[0].legend(handlelength=2, scatterpoints=1, numpoints=1,frameon=False,ncol=3,loc='lower center', bbox_to_anchor=(0.5, 1.0))
plt.savefig(sys.argv[1][:-4]+'.action.pdf',bbox_inches='tight')
plt.clf()

f,a=plt.subplots(1,1,figsize=(3,4))
plt.subplots_adjust(hspace=0.3)
a.plot(R,results.T[18],lw=0.5,color='k')
a.set_ylabel(r'$z$')
a.set_xlabel(r'$R$')
a.set_aspect('equal')
plt.savefig(sys.argv[1][:-4]+'.orbit.pdf',bbox_inches='tight')
