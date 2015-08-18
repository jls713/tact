import numpy as np
import matplotlib.pyplot as plt

Conv = 0.97775

g = np.genfromtxt("sag_wfreq_joined.dat",skip_footer=1)

f,a = plt.subplots(1,3,figsize=[10.,3.32])
for i in a:
	i.set_xlim(0.,2.)
	i.set_ylim(0.,2.)

Color = g.T[-2]

a[0].scatter(g.T[3]/np.pi,g.T[4]/np.pi,c = Color, vmin=min(Color), vmax=max(Color),lw=0, alpha=0.05)
a[0].set_xlabel(r'$\theta_1/\pi$')
a[0].set_ylabel(r'$\theta_2/\pi$')
a[1].scatter(g.T[3]/np.pi,g.T[5]/np.pi,c = Color, vmin=min(Color), vmax=max(Color),lw=0, alpha=0.05)
a[1].set_xlabel(r'$\theta_1/\pi$')
a[1].set_ylabel(r'$\theta_3/\pi$')
a[2].scatter(g.T[4]/np.pi,g.T[5]/np.pi,c = Color, vmin=min(Color), vmax=max(Color),lw=0, alpha=0.05)
a[2].set_xlabel(r'$\theta_2/\pi$')
a[2].set_ylabel(r'$\theta_3/\pi$')

plt.savefig("sag_angles.pdf")
plt.clf()

f,a = plt.subplots(1,3,figsize=[10.,3.32])
plt.subplots_adjust(wspace=0.3)
a[0].set_xlim(0.,4000.)
a[0].set_ylim(0.,2000.)
a[1].set_xlim(0.,4000.)
a[1].set_ylim(2000.,9000.)
a[2].set_xlim(0.,2000.)
a[2].set_ylim(2000.,9000.)

a[0].scatter(g.T[0],g.T[1],c = Color, vmin=min(Color), vmax=max(Color),lw=0, alpha=0.05)
a[0].set_xlabel(r'$J_1/\,{\rm kpc}\,{\rm kms}^{-1}$')
a[0].set_ylabel(r'$J_2/\,{\rm kpc}\,{\rm kms}^{-1}$')
a[1].scatter(g.T[0],g.T[2],c = Color, vmin=min(Color), vmax=max(Color),lw=0, alpha=0.05)
a[1].set_xlabel(r'$J_1/\,{\rm kpc}\,{\rm kms}^{-1}$')
a[1].set_ylabel(r'$J_3/\,{\rm kpc}\,{\rm kms}^{-1}$')
a[2].scatter(g.T[1],g.T[2],c = Color, vmin=min(Color), vmax=max(Color),lw=0, alpha=0.05)
a[2].set_xlabel(r'$J_2/\,{\rm kpc}\,{\rm kms}^{-1}$')
a[2].set_ylabel(r'$J_3/\,{\rm kpc}\,{\rm kms}^{-1}$')

f,a = plt.subplots(1,3,figsize=[10.,3.32])
plt.subplots_adjust(wspace=0.3)
a[0].set_xlim(2.,12.)
a[0].set_ylim(2.,8.)
a[1].set_xlim(2.,12.)
a[1].set_ylim(2.,9.)
a[2].set_xlim(2.,8.)
a[2].set_ylim(2.,9.)

a[0].scatter(g.T[6]/Conv,g.T[7]/Conv,c = Color, vmin=min(Color), vmax=max(Color),lw=0, alpha=0.05)
a[0].set_xlabel(r'$\Omega_1/\,{\rm Gyr}^{-1}$')
a[0].set_ylabel(r'$\Omega_2/\,{\rm Gyr}^{-1}$')
a[1].scatter(g.T[6]/Conv,g.T[8]/Conv,c = Color, vmin=min(Color), vmax=max(Color),lw=0, alpha=0.05)
a[1].set_xlabel(r'$\Omega_1/\,{\rm Gyr}^{-1}$')
a[1].set_ylabel(r'$\Omega_3/\,{\rm Gyr}^{-1}$')
a[2].scatter(g.T[7]/Conv,g.T[8]/Conv,c = Color, vmin=min(Color), vmax=max(Color),lw=0, alpha=0.05)
a[2].set_xlabel(r'$\Omega_2/\,{\rm Gyr}^{-1}$')
a[2].set_ylabel(r'$\Omega_3/\,{\rm Gyr}^{-1}$')

plt.savefig("sag_freq.pdf")
