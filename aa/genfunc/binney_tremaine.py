import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm
from genfunc_3d import find_actions, check_angle_solution
import test_potentials as pot

class logarithmic(object):
    def __init__(self):
        self.qz2 = 0.7**2
        self.qy2 = 0.9**2

    def pot(self,x):
        return 0.5*np.log(x[0]**2+x[1]**2/self.qy2+x[2]**2/self.qz2+0.1)

    def H(self,x):
    	return 0.5*(x[3]**2+x[4]**2+x[5]**2)+self.pot(x[:3])

    def tot_force(self,x,y,z):
        p = (x**2+y**2/self.qy2+z**2/self.qz2+0.1)
        return -np.array([x,y/self.qy2,z/self.qz2])/p

    def find_x(self,y,z):
    	return np.sqrt(np.exp(1.)-(y**2/self.qy2+z**2/self.qz2+0.1))

from numpy.fft import fft

def second_diff(x):
	sd = np.ones(len(x)-1,dtype=complex)
	sd[1:] = x[2:]+x[:-2]-2.*x[1:-1]
	sd[0] = 2.*(x[1]-x[0])
	return sd

def find_freqs(timeseries,results):
	freq = np.ones(3)
	n = len(timeseries)
	frq = np.arange(n)/timeseries[-1]
	frq = frq[range(n/2)]
	for i in range(3):
		X = fft(results.T[i])[range(n/2)]
		X[0] = np.real(X[0])
		X[-1] = np.real(X[-1])
		sd = second_diff(X)
		peak = np.argmax(np.abs(X))
		while(peak==0):
			print "H"
			X[0]+=sd[0]*np.pi/n
			sd = second_diff(X)
			peak = np.argmax(np.abs(X))
		sm,s,sp = sd[peak-1:peak+2]
		xm = -np.real((s+2.*sm)/(s-sm))
		xp = np.real((s+2.*sp)/(s-sp))
		x = .5*(xp+xm)
		freq[i]=(peak-x)*2.*np.pi/timeseries[-1]
	return freq


def plot_spectrum(timeseries,x,axis=None,xlabel="",ylabel=""):
	n = len(timeseries) # length of the signal
 	k = np.arange(n)
 	T = timeseries[-1]
 	frq = k/T # two sides frequency range
 	frq = frq[range(n/2)] # one side frequency range
 	X = fft(x)/n # fft computing and normalization
 	X = X[range(n/2)]
 	print(np.fft.fftfreq(timeseries.shape[-1]))
 	print(frq[1])
	if(axis==None):
		plt.plot(frq,np.abs(X))
		plt.show()
	else:
		axis.plot(frq,np.abs(X),'k')
		axis.set_xlabel(xlabel)
		axis.set_ylabel(ylabel)


def plots(timeseries,results):
	fig,axi = plt.subplots(2,3,figsize=[6.64,6.64/1.6])
	plt.subplots_adjust(wspace=0.4,hspace=0.25)
	axi[0][0].plot(results.T[0],results.T[1],'k')
	axi[0][0].set_xlabel(r'$x$')
	axi[0][0].set_ylabel(r'$y$')
	axi[0][1].plot(results.T[0],results.T[2],'k')
	axi[0][1].set_xlabel(r'$x$')
	axi[0][1].set_ylabel(r'$z$')
	axi[0][2].plot(results.T[1],results.T[2],'k')
	axi[0][2].set_xlabel(r'$y$')
	axi[0][2].set_ylabel(r'$z$')
	plot_spectrum(timeseries,results.T[0],axi[1][0],xlabel=r'Freq/$2\pi$',ylabel=r'$|$FT($x$)$|$')
	for i in range(3):
		axi[1][i].set_xlim(0.,0.4)
		axi[1][i].xaxis.set_major_locator(MaxNLocator(5))
		axi[0][i].xaxis.set_major_locator(MaxNLocator(5))
	plot_spectrum(timeseries,results.T[1],axi[1][1],xlabel=r'Freq/$2\pi$',ylabel=r'$|$FT($y$)$|$')
	plot_spectrum(timeseries,results.T[2],axi[1][2],xlabel=r'Freq/$2\pi$',ylabel=r'$|$FT($z$)$|$')
	plt.savefig("../../Documents/thesis/genfunc_aa/nonchaotic_spectrum.pdf",bbox_inches='tight')
	plt.clf()

from genfunc_3d import eval_mean_error_functions, check_angle_solution as ctas

def print_freqs():
	log = logarithmic()
	timeseries = np.linspace(0.,500.,2048)
	ntheta,nphi = 30,30
	for i in np.arange(4.,ntheta+1)/(ntheta+1):
		for j in .5*np.pi*np.arange(15.,nphi+1)/(nphi+1):
			# if(i>0.9 and j>np.pi*0.45):
			# if(j>np.pi/4.):
			# if(i>0.59 and j>np.pi/2.*0.9):
			print np.arccos(i),j
			st,ct = np.sin(np.arccos(i)),i
			sp,cp = np.sin(j),np.cos(j)
			combo = cp*cp*st*st+sp*sp*st*st/log.qy2+ct*ct/log.qz2
			r = np.sqrt((np.exp(1.)-0.1)/combo)
			initial = np.array([r*cp*st,r*sp*st,r*ct,0.0001,0.0001,0.0001])
			# initial = np.array([0.111937987197,0.0104758765442,1.12993449025,0.0001,0.0001,0.0001])
			results = odeint(pot.orbit_derivs2,initial,timeseries,args=(log,),rtol=1e-5,atol=1e-5)
			# print(log.H(initial),log.H(results[-1]))
			plots(timeseries,results)
			# freq = find_freqs(timeseries,results)
			L = find_actions(results, timeseries, N_matrix = 4, ifloop=True,ifprint = False)
			# if(L==None):
				# break
			(act,ang,n_vec,toy_aa,para),loop = L
			E = eval_mean_error_functions(act,ang,n_vec,toy_aa,timeseries,withplot=False)/np.std(timeseries)
			# ctas(ang,n_vec,toy_aa,timeseries)
			# print freq[0],freq[1],freq[2]
			print ang[3],ang[4],ang[5],initial[0],initial[1],initial[2],act[0],act[1],act[2] #,E[0],E[1],E[2],E[3],E[4],E[5] #,freq[0],freq[1],freq[2]
			exit()

def plot_freqs():
	plt.figure(figsize=[3.32,3.6])
	g = np.genfromtxt("with_fft")
	omega21 = map(lambda i,j: np.abs(j)/i*2. if i>1.2 else j/i,g.T[0],g.T[1])
	omega31 = map(lambda i,j: j/i*2. if 
		i>1.2 else j/i,g.T[0],g.T[2])
	# G = np.sqrt(g.T[12]**2) #+g.T[13]**2+g.T[14]**2)
	# C = G/np.max(G)	
	plt.scatter(omega21,omega31,c='k',marker ='.',s=1.,edgecolors="none")
	plt.xlabel(r'$\Omega_2/\Omega_1$')
	plt.ylabel(r'$\Omega_3/\Omega_1$')
	plt.ylim(1.1,2.1)
	plt.savefig('binney_tremaine_fig345.pdf')

print_freqs()
# plot_freqs()
