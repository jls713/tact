# Solving the series of linear equations for true action
# and generating function Fourier components

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import fmin
from scipy.integrate import quad


def sho_H(x,omega):
    """ Simple harmonic oscillator Hamiltonian = 0.5 * omega**2 * x**2"""
    return 0.5*(x[1]**2+(omega*x[0])**2)


def sho_derivs(x,t,omega):
    """ Derivatives of sho pot for orbit integration """
    return np.array([x[1],-omega**2*x[0]])


def other_H(x):
    """ Quartic potential Hamiltonian """
    return 0.5*x[1]**2+0.25*x[0]**4


def otherderivs(x,t):
    """ Derivatives of quartic potential for orbit integration """
    return np.array([x[1],-x[0]**3])


def minfnc(omega,xsamples):
    Htrue=other_H(xsamples[0])
    return np.sum([(sho_H(i,omega[0])-Htrue)**2 for i in xsamples])


def minfnc_deriv(omega,xsamples):
    Htrue=other_H(xsamples[0])
    return np.sum([i[0]**2*(sho_H(i,omega[0])-Htrue) for i in xsamples])


def findbestomega(xsamples):
    """ Minimize sum of square differences between H and H_sho"""
    return fmin(minfnc,np.array([1.]), args=(xsamples,))


def trueaction(x):
    """ Find true action for quartic potential \Phi = 0.25*x**4 """
    En = other_H(x)
    xlim=(4.*En)**0.25
    return 2.*quad(lambda y:np.sqrt(2.*En-0.5*y**4),0.,xlim)[0]/np.pi


def truefreq(x):
    """ Find true freq. for quartic potential \Phi = 0.25*x**4 """
    En = other_H(x)
    xlim=(4.*En)**0.25
    return np.pi/quad(lambda y:2./np.sqrt(2.*En-0.5*y**4),0.,xlim)[0]


def angact(x,omega):
    """ Calculate angle and action variable in sho potential with
    parameter omega """
    action = (x[1]**2+(omega*x[0])**2)/(2.*omega)
    angle = np.arctan(-x[1]/omega/x[0])
    if(x[0]<0):
        angle+=np.pi
    return np.array([action,angle % (2.*np.pi)])

from numpy.linalg import solve


def solver(AA, N_coeffs):
    n = N_coeffs+1
    b = np.zeros(shape=(n, ))
    a = np.zeros(shape=(n,n))
    a[0,0]=len(AA)

    for i in AA:
        for k in range(n):
            for l in range(k,n):
                if k==0:
                    a[k,l]+=2.*l*(np.cos(l*i[1]))
                else:
                    a[k,l]+=4.*k*l*(np.cos(k*i[1])*np.cos(l*i[1]))
        b[0]+=i[0]
        for j in range(1,n):
            b[j]+=2.*j*i[0]*(np.cos(j*i[1]))

    a = (a+a.T)
    for i in range(n):
        a[i,i]/=2.

    return solve(a,b)

from itertools import izip


def angle_solver(AA, timeseries, N_coeffs):
    n=N_coeffs+2
    b = np.zeros(shape=(n, ))
    a = np.zeros(shape=(n,n))
    a[0,0]=len(AA)
    a[0,1]=np.sum(timeseries)
    a[1,1]=np.sum(timeseries*timeseries)
    for i,j in izip(AA,timeseries):
        for k in range(2,n):
            for l in range(0,k+1):
                if l==0:
                    a[k,l]+=-2.*np.sin((k-1)*i[1])
                elif l==1:
                    a[k,l]+=-2.*j*np.sin((k-1)*i[1])
                else:
                    a[k,l]+=4.*np.sin((k-1)*i[1])*np.sin((l-1)*i[1])
        b[0]+=i[1]
        b[1]+=j*i[1]
        for k in range(2,n):
            b[k]+=-2.*i[1]*np.sin((k-1)*i[1])

    a = (a+a.T)
    for i in range(n):
        a[i,i]/=2.

    return solve(a,b)


def potential_plot(omega):
    x = np.linspace(-1,1,100)
    plt.plot(x,0.5*omega*x*x,label='toy')
    plt.plot(x,0.25*x*x*x*x,label='target')
    plt.legend()
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\Phi(x)$')
    plt.show()


def plotter(omega,N_coeffs,N_samples):

    timeseries=np.linspace(0.,7.,N_samples)
    results = odeint(otherderivs,np.array([1.,0.]),timeseries)
    if(omega==0):
        omega = findbestomega(results)[0]
    print omega
    potential_plot(omega)
    AA = np.array([angact(i,omega) for i in results])
    Energy = np.array([sho_H(i,omega) for i in results])
    Energy_exact = np.array([other_H(i) for i in results])
    actions=np.array([trueaction(i) for i in results])

    print 'True action = ',actions[0]
    print 'Average action = ',np.mean(AA.T[0])
    solved = solver(AA,N_coeffs)
    print 'Estimated action = ',solved[0]

    f,a=plt.subplots(6,1)
    a[0].plot(timeseries,results.T[0])
    a[0].set_ylabel(r'$x$')
    a[1].plot(timeseries,results.T[1])
    a[1].set_ylabel(r'$v$')
    a[2].plot(timeseries,Energy,label='toy energy')
    a[2].plot(timeseries,Energy_exact, label='exact energy')
    a[2].set_ylabel(r'$E$')
    a[2].legend()
    a[3].plot(timeseries,AA.T[0],label='toy action')
    a[3].plot(timeseries,actions,label='true action')
    a[3].plot(timeseries,solved[0]*np.ones(len(timeseries)),label='est. action')
    a[3].legend()
    a[3].set_ylabel(r'$J$')
    a[4].plot(timeseries,AA.T[1],label='toy angle')
    a[4].set_ylabel(r'$\theta$')
    a[4].set_xlabel(r'$t$')
    a[4].legend()
    g = np.arange(1,N_coeffs+1)
    a[5].plot(g,np.log10(np.abs(solved[1:])),'.',label=r'$\log_{10}\Delta J=$'+str(np.log10(np.abs(solved[0]-actions[0]))))
    a[5].set_ylabel(r'$\log_{10}|S_n|$')
    a[5].set_xlabel(r'$n$')
    a[5].legend()
    plt.show()


def delta_act(omega,N_coeffs,N_samples):

    timeseries=np.linspace(0.,7.5,N_samples)
    results = odeint(otherderivs,np.array([1.,0.]),timeseries)
    if(omega==0):
        omega = findbestomega(results)[0]
    AA = np.array([angact(i,omega) for i in results])
    avact=np.mean(AA.T[0])
    action=trueaction(results[0])
    solved = solver(AA,N_coeffs)
    return np.abs(solved[0]-action), np.abs(avact-action)


def delta_act_plots():
    n_coeff = np.arange(1,50)
    D = np.array([delta_act(0,i,50) for i in n_coeff])
    D2 = np.array([delta_act(0,i,100) for i in n_coeff])
    plt.plot(n_coeff,np.log10(D),'.',label='50 time samples, best omega')
    plt.plot(n_coeff,np.log10(D2),'.',label='100 time samples, best omega')
    plt.xlabel(r'# of S_n')
    plt.ylabel(r'$\log_{10}|J-J_{\rm true}|$')
    plt.legend()
    plt.show()

fig,a = plt.subplots(2,1)
x = np.linspace(-1,1,100)
for i in np.linspace(0.3,3.,50):
    Da,Aa = delta_act(i,25,50)
    a[0].plot(x,0.5*i*i*x*x,c='k',alpha=0.2)
    a[1].plot(i,np.log10(Da),'.',c='k')
    a[1].plot(i,np.log10(Aa),'.',c='b')

a[0].plot(x,0.25*x*x*x*x,c='r')
a[1].axvline(0.63544921875,c='r')
a[0].set_xlabel(r'$x$')
a[0].set_ylabel(r'$\Phi(x)$')
a[1].set_xlabel(r'$\omega$')
a[1].set_ylabel(r'$\log_{10}\Delta J$')
plt.show()
# timeseries=np.linspace(0.,7.,50)
# results = odeint(otherderivs,np.array([1.,0.]),timeseries)
# omega = findbestomega(results)[0]
# AA = np.array([angact(i,omega) for i in results])
# #fig,a = plt.subplots(2,1)
# omexact = truefreq(results[0])
# a=np.array([angle_solver(AA,timeseries,i)[1] for i in range(50)])
# plt.plot(range(50),np.log10(np.abs(a-omexact)),'.',c='k')
# plt.xlabel(r'# of $S_n$')
# plt.ylabel(r'$\log_{10}|\Omega-\Omega_{\mathrm{true}}|$')
# # angs = SOL[0]+SOL[1]*timeseries
# # snarray=np.arange(1,21)
# # a[0].plot(timeseries,angs,c='k',label='target')
# # a[0].plot(timeseries,AA.T[1],c='r',label='toy')
# # a[0].legend(loc=2)
# # a[0].set_xlabel(r'$t$')
# # a[0].set_ylabel(r'$\theta$')
# # a[1].plot(snarray,np.log10(np.abs(SOL[2:])),'.',c='k')
# # a[1].set_xlabel(r'$n$')
# # a[1].set_ylabel(r'$\log_{10} |dS_n/dJ\'|$')
# plt.show()
# # print SOL[0]%(2.*np.pi)
# #plotter(0,25,50)
