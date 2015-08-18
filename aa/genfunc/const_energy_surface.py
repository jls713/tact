import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm
from genfunc_3d import find_actions, check_angle_solution
from scipy.optimize import brentq

import test_potentials as pot
import solver

# For plot used YMax = 18, theta_min = 0.05, 10 by 10 , t_max=10


LM = pot.LMPot()

def rotate_coords(X):
	xy = LM.invcoordrot(X[0],X[1])
	vxvy = LM.invcoordrot(X[3],X[4])
	return np.array([xy[0],xy[1],X[2],vxvy[0],vxvy[1],X[5]])

def unrotate_coords(X):
	xy = LM.coordrot(X[0],X[1])
	vxvy = LM.coordrot(X[3],X[4])
	return np.array([xy[0],xy[1],X[2],vxvy[0],vxvy[1],X[5]])

zlittle = 0.05
YMax = 18.
xymax = LM.coordrot(0.,YMax)
y = np.linspace(0.2,YMax-0.2,50)
E_cons = LM.tot_pot(xymax[0],xymax[1],zlittle)
# print(E_cons)
timeseries=np.linspace(0.,10.,1000)
# timeseries_long=np.linspace(0.,50.,500)

from solver import check_each_direction as ced
from solver import unroll_angles as ua
from genfunc_3d import assess_angmom, eval_mean_error_functions

for i in y:
	# i = 8.11111111111
	# i = 1.27755102041
	xy = LM.coordrot(0.,i)
	Phi = LM.tot_pot(xy[0],xy[1],zlittle)
	pylittle = 0.01*np.sqrt(2.*(E_cons-Phi))
	P = np.sqrt(2.*(E_cons-Phi)-pylittle**2)
	for NN,theta in enumerate(np.linspace(0.05,np.pi/2.-0.05,50)):
		px = P*np.cos(theta)
		pz = P*np.sin(theta)
		# px,pz =  116.216249287,200.519902867 
		# px,pz = 156.389296517, 471.688963618
		pxpy = LM.coordrot(px,pylittle)
		initial = np.array([xy[0],xy[1],zlittle,pxpy[0],pxpy[1],pz])
		results = np.array([rotate_coords(p) for p in odeint(pot.orbit_derivs2,initial,timeseries,args=(LM,),rtol=1e-10,atol=1e-10)])
		# results_long = np.array([rotate_coords(p) for p in odeint(pot.orbit_derivs2,initial,timeseries_long,args=(LM,),rtol=1e-5,atol=1e-5)])
		# loop = assess_angmom(results_long)	
		# LL = np.any(loop>0)
		# BOX = not LL
		L = find_actions(results, timeseries, N_matrix = 6, ifloop=True,ifprint = False)

		if(L==None):
			print LM.H(initial),i,px,pz,loop
			continue

		(act,ang,n_vec,toy_aa,para),loop = L
		# print(eval_mean_error_functions(act,ang,n_vec,toy_aa,timeseries,withplot=True))
		checks,maxgap = ced(n_vec,ua(toy_aa.T[3:].T,np.ones(3)))
		TMAX = 1000
		counter = 0
		while(len(checks)>0 and counter<10):
			# print(checks)
			if(np.any(maxgap<0.)):
				TMAX=TMAX*2.
				results = np.array([rotate_coords(p) for p in odeint(pot.orbit_derivs2,initial,np.linspace(0.,10.,TMAX),args=(LM,),rtol=1e-10,atol=1e-10)])
				L = find_actions(results, np.linspace(0.,10.,TMAX), N_matrix = 6, ifloop=True,ifprint = False)
				print("Finer sampling")
			else:
				results2 = np.array([rotate_coords(p) for p in odeint(pot.orbit_derivs2,unrotate_coords(results[-1]),timeseries,args=(LM,),rtol=1e-10,atol=1e-10)])
				results = np.vstack((results,results2))
				print("Longer integration")
		# 	print(results.T[-1])
				L = find_actions(results, np.linspace(0.,len(results)/100.,len(results)), N_matrix = 6, ifloop=True,ifprint = False)
				# print(len(results))
			(act,ang,n_vec,toy_aa,para),loop = L
			print(len(results),loop)
			checks,maxgap = ced(n_vec,ua(toy_aa.T[3:].T,np.ones(3)))
			counter=counter+1
		# if(counter==10):
			# continue
		# print(ang)
		minfreq = np.min(np.abs(ang[3:6]))
		if(8.*2.*np.pi/minfreq>10.):
			print(i,px,pz,"Fewer than 8 fundamental periods")
		T = np.linspace(0.,len(results)/100.,len(results))
		errors = eval_mean_error_functions(act,ang,n_vec,toy_aa,T,withplot=False)/np.sqrt(len(T))
		print LM.H(initial),i,px,pz,act[0],act[1],act[2],' '.join(map(str, loop)),ang[3],ang[4],ang[5], errors[0],errors[1], LM.H(unrotate_coords(results[-1]))
		



		# check_angle_solution(ang,n_vec,toy_aa,timeseries)
		# plt.plot(toy_aa.T[3],toy_aa.T[4],'.')
		# plt.show()
		# plt.plot(toy_aa.T[3],toy_aa.T[5],'.')
		# plt.show()
		# plt.plot(results.T[0],results.T[1],'.');plt.show()
		# F = pot.orbit_integrate(initial,80.,LM)
		# results = np.array([rotate_coords(p) for p in F[0]])
		# print(len(results))
		# plt.plot(results.T[0],results.T[1])
		# plt.show()
		# timeseries = F[1]
		# plt.plot(np.sqrt(results.T[0]**2+results.T[1]**2),results.T[2])
		# plt.show()

# timeseries=np.linspace(0.,15.,500)
# results = odeint(pot.orbit_derivs,np.array([1.21,1.,0.6,200.,200.,200.]),timeseries,args=(LM,))
# print LM.H(results[0])
# plt.plot(results.T[0],results.T[1])
# plt.show()
