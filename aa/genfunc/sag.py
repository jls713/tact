import numpy as np
from scipy.integrate import odeint
import test_potentials as pot
from genfunc_3d import find_actions
import matplotlib.pyplot as plt
from solver import check_each_direction as ced, unroll_angles as ua

data = np.genfromtxt("SgrTriax_DYN.dat")
datae = np.vstack((data.T[6],data.T[7],data.T[8],data.T[15],data.T[16],data.T[17])).T[1:]

timeseries = np.linspace(0.,10.,1000)

LM = pot.LMPot()

def rotate_coords(X):
	xy = LM.invcoordrot(X[0],X[1])
	vxvy = LM.invcoordrot(X[3],X[4])
	return np.array([xy[0],xy[1],X[2],vxvy[0],vxvy[1],X[5]])


def unrotate_coords(X):
	xy = LM.coordrot(X[0],X[1])
	vxvy = LM.coordrot(X[3],X[4])
	return np.array([xy[0],xy[1],X[2],vxvy[0],vxvy[1],X[5]])


for i in datae:
	results = np.array([rotate_coords(p) for p in odeint(pot.orbit_derivs2,i,timeseries,args=(LM,),rtol=1e-10,atol=1e-10)])
	L = find_actions(results,timeseries,N_matrix=6,ifprint=False)
	(act,ang,n_vec,toy_aa,para) = L
	checks,maxgap = ced(n_vec,ua(toy_aa.T[3:].T,np.ones(3)))
	TMAX = 1000
	counter = 0
	while(len(checks)>0 and counter<10):
		# print(checks)
		if(np.any(maxgap<0.)):
			TMAX=TMAX*2.
			results = np.array([rotate_coords(p) for p in odeint(pot.orbit_derivs2,initial,np.linspace(0.,10.,TMAX),args=(LM,),rtol=1e-10,atol=1e-10)])
			L = find_actions(results, np.linspace(0.,10.,TMAX), N_matrix = 6, ifprint = False)
		else:
			results2 = np.array([rotate_coords(p) for p in odeint(pot.orbit_derivs2,unrotate_coords(results[-1]),timeseries,args=(LM,),rtol=1e-10,atol=1e-10)])
			results = np.vstack((results,results2))
			L = find_actions(results, np.linspace(0.,len(results)/100.,len(results)), N_matrix = 6, ifprint = False)
		(act,ang,n_vec,toy_aa,para) = L
		checks,maxgap = ced(n_vec,ua(toy_aa.T[3:].T,np.ones(3)))
		counter+=1
	(act,ang,n_vec,toy_aa,para) = L
	print act[0],act[1],act[2],ang[0],ang[1],ang[2],ang[3],ang[4],ang[5],i[-2],i[-1]