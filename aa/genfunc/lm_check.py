import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm
from genfunc_3d import find_actions

# in units kpc, km/s and 10^11 M_solar
Grav = 430091.5694
Conv = 0.97775

import test_potentials as pot
import solver

x = 0.
z = 0.
py = 0.

LM = pot.LMPot()

def rotate_coords(X):
	xy = LM.invcoordrot(X[0],X[1])
	vxvy = LM.invcoordrot(X[3],X[4])
	return np.array([xy[0],xy[1],X[2],vxvy[0],vxvy[1],X[5]])

zlittle = 0.1234234
YMax = 18.
y=14.8
xymax = LM.coordrot(0.,YMax)
LowPer = 2.*np.pi/18.229468031
E_cons = LM.tot_pot(xymax[0],xymax[1],zlittle)
timeseries=np.linspace(0.,LowPer*16.,1000)
# timeseries=np.linspace(0.,LowPer*8.,2000)
xy = LM.coordrot(0.,y)
Phi = LM.tot_pot(xy[0],xy[1],zlittle)
pz = 0.32529555234*np.sqrt(2.*(E_cons-Phi))
pylittle = 0.001*np.sqrt(2.*(E_cons-Phi))
px = 2.*(E_cons-0.5*(pz*pz+pylittle**2)-Phi)
px = -np.sqrt(px)
pxpy = LM.coordrot(px,pylittle)
initial = np.array([xy[0],xy[1],zlittle,pxpy[0],pxpy[1],pz])
print initial
# initial = np.array([14.7,1.8,0.1,16.0,-128.9,44.7])
print LM.H(initial)
results_LM = odeint(pot.orbit_derivs2,initial,timeseries,args=(LM,))
results = np.array([rotate_coords(p) for p in results_LM])

Roll = 500

f,a = plt.subplots(2,3,figsize=[3.32,5.5])
a[0,0] = plt.subplot2grid((3,2), (0, 0))
a[1,0] = plt.subplot2grid((3,2), (0, 1))
a[0,1] = plt.subplot2grid((3,2), (1, 0))
a[1,1] = plt.subplot2grid((3,2), (1, 1))
a[0,2] = plt.subplot2grid((3,2), (2, 0))
a[1,2] = plt.subplot2grid((3,2), (2, 1))
plt.subplots_adjust(wspace=0.6,hspace=0.45)

a[0,0].set_xlabel(r'$x/{\rm kpc}$')
a[1,0].set_xlabel(r'$x/{\rm kpc}$')
a[0,0].set_ylabel(r'$y/{\rm kpc}$')
a[1,0].set_ylabel(r'$z/{\rm kpc}$')
a[0,1].set_xlabel(r'$t/{\rm Gyr}$')
a[0,1].set_ylabel(r"$(J_1'-\langle J_1'\rangle )/{\rm kpc\,km\,s}^{-1}$")
a[1,1].set_xlabel(r'$t/{\rm Gyr}$')
a[1,1].set_ylabel(r"$(\Omega_1'-\langle\Omega_1'\rangle )/1000\,{\rm Gyr}^{-1}$")
a[0,2].set_xlabel(r'$\theta_1/\pi$')
a[0,2].set_ylabel(r'$\theta_2/\pi$')
a[1,2].set_xlabel(r'$\theta_1/\pi$')
a[1,2].set_ylabel(r'$\theta_3/\pi$')

times = np.array([])
angles = np.array([])
toyangles = np.array([])
actions = np.array([])
freqs= np.array([])
(act,ang,n_vec,toy_aa,para),loop = find_actions(results, timeseries, N_matrix = 6, ifloop=True,ifprint = False)
size = len(ang[6:])/3
AA = [np.array([np.sum(ang[6+i*size:6+(i+1)*size]*np.sin(np.sum(n_vec*K,axis=1))) for K in toy_aa.T[3:].T]) for i in range(3)]
# a[0,2].plot((toy_aa.T[3]+2.*AA[0]) % (2.*np.pi) / np.pi,(toy_aa.T[4]+2.*AA[1]) % (2.*np.pi)/ np.pi,'.',markersize=3)
# a[1,2].plot((toy_aa.T[3]+2.*AA[0]) % (2.*np.pi) / np.pi,(toy_aa.T[5]+2.*AA[2]) % (2.*np.pi) / np.pi,'.',markersize=3)
# plt.savefig('LM_test.pdf',bbox_inches='tight')

# exit(0)
toyangles = toy_aa.T[3:6].T/np.pi
MaxRoll = 0.5*Roll
from solver import check_each_direction as ced
from solver import unroll_angles as ua

for i in np.arange(0,MaxRoll,Roll/100.):
	(act,ang,n_vec,toy_aa,para),loop = find_actions(results[i:i+Roll], timeseries[:Roll], N_matrix = 6, ifloop=True,ifprint = False)
	checks,maxgap = ced(n_vec,ua(toy_aa.T[3:].T,np.ones(3)))
	print(checks)
	print len(results[i:i+Roll]),act[0],act[1],act[2],ang[3],ang[4],ang[5],ang[0] % (2.*np.pi),ang[1] % (2.*np.pi),ang[2] % (2.*np.pi),' '.join(map(str, loop))
	if(len(angles)==0):
		angles = ang[:3]%(2.*np.pi)/np.pi
		freqs = ang[3:6]
		actions = act[:3]
		times = timeseries[i]
	else:
		angles = np.vstack((angles,ang[:3]%(2.*np.pi)/np.pi))
		freqs = np.vstack((freqs,ang[3:6]))
		actions = np.vstack((actions,act[:3]))
		times = np.vstack((times,timeseries[i]))
		
freqs = freqs/Conv
print(np.mean(actions.T[0]),np.mean(actions.T[1]),np.mean(actions.T[2]))
print(np.mean(freqs.T[0]),np.mean(freqs.T[1]),np.mean(freqs.T[2]))
print(np.std(actions.T[0]),np.std(actions.T[1]),np.std(actions.T[2]))
print(np.std(freqs.T[0]),np.std(freqs.T[1]),np.std(freqs.T[2]))
a[0,0].plot(results_LM.T[0],results_LM.T[1],'k')
a[1,0].plot(results_LM.T[0],results_LM.T[2],'k')
a[0,1].plot(Conv*times,actions.T[0]-np.mean(actions.T[0]),'k.',markersize=3)
a[1,1].plot(Conv*times,(freqs.T[0]/Conv-np.mean(freqs.T[0])/Conv)*1000.,'k.',markersize=3)
a[0,2].plot(angles.T[0],angles.T[1],'k.',zorder=1,markersize=3)
a[1,2].plot(angles.T[0],angles.T[2],'k.',zorder=1,markersize=3)
# a[0,2].plot(toyangles.T[0][:MaxRoll],toyangles.T[1][:MaxRoll],'r.',linewidth=0.2,markersize=0.4,alpha=0.4,zorder=0)
a[0,2].plot((timeseries[:MaxRoll]*freqs[0][0]+angles[0][0]*np.pi)%(2.*np.pi)/np.pi,(timeseries[:MaxRoll]*freqs[0][1]+angles[0][1]*np.pi) % (2.*np.pi)/np.pi,'b.',linewidth=0.2,markersize=0.4,alpha=0.7,zorder=0)
# a[1,2].plot(toyangles.T[0][:MaxRoll],toyangles.T[2][:MaxRoll],'r.',linewidth=0.2,markersize=0.4,alpha=0.4,zorder=0)
a[1,2].plot((timeseries[:MaxRoll]*freqs[0][0]+angles[0][0]*np.pi)%(2.*np.pi)/np.pi,(timeseries[:MaxRoll]*freqs[0][2]+angles[0][2]*np.pi)%(2.*np.pi)/np.pi,'b.',linewidth=0.2,markersize=0.4,alpha=0.7,zorder=0)
a[0,1].xaxis.set_major_locator(MaxNLocator(4))
a[1,1].xaxis.set_major_locator(MaxNLocator(4))
a[0,0].xaxis.set_major_locator(MaxNLocator(5))
a[1,0].xaxis.set_major_locator(MaxNLocator(5))
plt.savefig('LM_test.pdf',bbox_inches='tight')
