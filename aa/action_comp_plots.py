### action_comp_plots.py
##  Produces the plots in Sanders & Binney (2016)
##  Fig 2. is produced by method_comparison which requires the result of running orbits.sh
##  Figs 3., 4., 5. and 6. are produced by time_plot which requires the result of running
##  mains/./many_tori.exe many_tori_output.dat
##  Fig 7. uses the output from orbits_converg.sh

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import lines
import sys
import seaborn as sns

# plt.rc('text', usetex=True)

kms2kpcGyr = 1000./977.775

def method_comparison(results_file,name,dontplot=[],rangess=None,Jranges=None):
	results = np.genfromtxt(results_file)
	labels = [r'Fudge v1','ItTC','O2GF','AvGF',r'Fudge v2','CAA','SAA',r'Fit']
	f = plt.figure(figsize=(6.64,3))

	a=[plt.subplot2grid((2,2),(i,0),colspan=1,rowspan=1) for i in range(2)]

	colors = [sns.color_palette()[2],sns.color_palette()[0],sns.color_palette()[1],sns.color_palette("Set2", 10)[3],sns.color_palette()[3],sns.color_palette()[5],sns.color_palette()[4],sns.color_palette("Set1",5,desat=.6)[4]]

	dashed=['-','--','--','--','-','-','-','-']

	times = results.T[0]

	plt.subplots_adjust(hspace=0.)
	JR = np.mean(results.T[5])
	Jz = np.mean(results.T[6])
	R = np.sqrt(results.T[17]*results.T[17]+results.T[18]*results.T[18])
	T=[0.,0.]
	NN=0
	for i in range(1,len(R)-2):
		if(R[i-1]>R[i] and R[i+1]>R[i]):
			if(T[0]==0.):
				T[0]=times[i]
			else:
				T[1]=times[i]
				NN+=1
			# for i in range(2):
	# 	a[i].set_xlim(0.,1000.)
	TR = (T[1]-T[0])/NN
	Mx = times[-1]/TR
	if(rangess):
		Mx=rangess
	spreads = np.zeros(shape=(9,2))
	spreads[0][0]=np.mean(results.T[5]) # pass actions as first argument
	spreads[0][1]=np.mean(results.T[6])
	for N,i in enumerate([5,6,0,4,7,1,3,2]):
		JRspread = np.std(results.T[i*2+1])
		Jzspread = np.std(results.T[i*2+2])
		## pass back std wrt genfunc method
		spreads[N+1]=np.array([np.sqrt(np.mean((results.T[i*2+1]-spreads[0][0])**2)),
		                       np.sqrt(np.mean((results.T[i*2+2]-spreads[0][1])**2))])
		n=i
		if(i>5):
			n=i-6
		if i in dontplot:
			continue
		l1, = a[0].plot((times-T[0])/TR,results.T[i*2+1],label=labels[i],lw=.5,ls=dashed[i],color=colors[i])
		l2, = a[1].plot((times-T[0])/TR,results.T[i*2+2],lw=0.5,ls=dashed[i],color=colors[i])
		if(i==1 or i==2 or i==3):
			l1.set_dashes((2,1))
			l2.set_dashes((2,1))
		x,y=np.array([[1.05*Mx+Mx*0.01*N,1.05*Mx+Mx*0.01*N],[JR-JRspread,JR+JRspread]])
		line = lines.Line2D(x, y, lw=1.,ls=dashed[i],color=colors[i])
		if(i==1 or i==2 or i==3):
			line.set_dashes((0.5,0.25))
		line.set_clip_on(False)
		a[0].add_line(line)
		x,y=np.array([[1.05*Mx+Mx*0.01*N,1.05*Mx+Mx*0.01*N],[Jz-Jzspread,Jz+Jzspread]])
		line = lines.Line2D(x, y, lw=1.,ls=dashed[i],color=colors[i])
		if(i==1 or i==2 or i==3):
			line.set_dashes((0.5,0.25))
		line.set_clip_on(False)
		a[1].add_line(line)
	if(Jranges):
		a[0].set_ylim(Jranges[0][0],Jranges[0][1])
		a[1].set_ylim(Jranges[1][0],Jranges[1][1])

	plt.setp(a[0].get_xticklabels(),visible=False)
	plt.setp(a[1].get_yticklabels()[-1],visible=False)

	a[0].set_xlabel(r'$t/T_R$')
	a[1].set_xlabel(r'$t/T_R$')
	a[0].set_ylabel(r'$J_R/\mathrm{kpc\,km\,s}^{-1}$')
	a[1].set_ylabel(r'$J_z/\mathrm{kpc\,km\,s}^{-1}$')

	leg=a[0].legend(handlelength=0.8, scatterpoints=1, numpoints=1,frameon=False,ncol=4,loc='lower left', bbox_to_anchor=(0.05, 1.0),fontsize=10)
	# set the linewidth of each legend object
	for legobj in leg.legendHandles:
	    legobj.set_linewidth(2.0)
	a[0].set_xlim(0.,Mx)
	a[1].set_xlim(0.,Mx)

	a = plt.subplot2grid((2,2),(0,1),colspan=1,rowspan=2)
	a.plot(R,results.T[19],lw=0.5,color='k')
	a.set_ylabel(r'$z/\mathrm{kpc}$')
	a.set_xlabel(r'$R/\mathrm{kpc}$')
	a.set_aspect('equal')
	a.yaxis.tick_right()
	a.yaxis.set_label_position("right")
	a.annotate(name, xy=(0.95,1.02), xycoords='axes fraction', fontsize=16,
	                horizontalalignment='right', verticalalignment='bottom')
	plt.savefig(results_file+'.action.pdf',bbox_inches='tight')
	return spreads

kappa = 45.7439
nu =  72.3666
## Piffl 2014
thin_sigmaR = 34.
thin_sigmaz = 25.
thick_sigmaR = 51.
thick_sigmaz = 49.

## Williams 2014
alpha0 = 0.84
alphainf=9.13
f = 0.59
delta = 0.01
v0=1.982
r0=2.16
J_b=3600.


def disc_dff(J,sigma,freq):
	return freq/sigma**2*np.ones(len(J))

def halo_dff(J,alpha0,alpha_inf,f,delta,v0,r0,rad=False):
	absJ=np.sqrt(np.sum(J**2))
	zeta = 2.*delta/(delta-2.)
	# J_b=1000.*np.power(delta/zeta,(1./delta))/np.sqrt(delta)*np.power(1.-0.5*delta*(1-1./np.e),(1./zeta))*v0*r0**(-delta/zeta)
	D0 = np.sqrt(2.*np.pi/np.e)
	curl = J.T[1]+J.T[2]+f*D0*J.T[0]
	lambd = zeta*(alpha0/delta-1.5)
	mu = zeta*(alphainf/delta-1.5)
	fac = 1.
	if(rad):
		fac = np.fabs(f*D0)
	return fac*np.fabs(-lambd/curl-(mu-lambd)*curl/(curl**2+J_b**2))

def find_JRJz_of_resonance(results):
	N=np.zeros(4)
	OmRes = np.array([.5,2./3.,1.,4./3.])
	for i in range(len(results)):
		if results[i][7]/results[i][9] < 1./2.:
			N[0]+=1
		if results[i][7]/results[i][9] < 2./3.:
			N[1]+=1
		if results[i][7]/results[i][9] < 1.:
			N[2]+=1
		if results[i][7]/results[i][9] < 4./3.:
			N[3]+=1
	return map(lambda i:(results[i[0]][3]-results[i[0]-1][3])/(results[i[0]][7]/results[i[0]][9]-results[i[0]-1][7]/results[i[0]-1][9])*(i[1]-results[i[0]-1][7]/results[i[0]-1][9])+results[i[0]-1][3], zip(N,OmRes))


def draw_res_lines(ax,res_line,yval):
	for i in range(len(ax)):
		for j in range(4):
			l = ax[i].axvline(res_line[j],color='k',ls='dashed',alpha=0.5)
			l.set_dashes((2,1))
	res_labels = ['1:2','2:3','1:1','4:3']
	# yval = ax[-1].get_ylim()[0]*2.
	for n,r in enumerate(res_labels):
		if(n==3):
			ax[-1].annotate(r,xy=(res_line[n]*0.95,yval), horizontalalignment='right', verticalalignment='bottom',fontsize=10)
		else:
			ax[-1].annotate(r,xy=(res_line[n]*1.05,yval), horizontalalignment='left', verticalalignment='bottom',fontsize=10)

def time_plot(results_file,with_res=False,with_boxes=False):
	results = np.genfromtxt(results_file)
	results = results[results.T[3]<0.6]
	res_line=0.
	if(with_res):
		res_line = find_JRJz_of_resonance(results)
	labels = [r'Fudge v1','ItTorus','O2GF','AvGF',r'Fudge v2','CAA','SAA',r'Fit']
	align = ['center','center','top','top','bottom','center','bottom','center']
	fig,a = plt.subplots(2,1,figsize=(4.,8.))

	plt.subplots_adjust(hspace=0.)
	colors = [sns.color_palette()[2],sns.color_palette()[0],sns.color_palette()[1],sns.color_palette("Set2", 10)[3],sns.color_palette()[3],sns.color_palette()[5],sns.color_palette()[4],sns.color_palette("Set1",5,desat=.6)[4]]

	dashed=['-','--','--','--','-','-','-','-']

	MaxX=0.6

	JRJz=results.T[3]
	JR = results.T[0]
	Jz = results.T[2]
	Jp = results.T[1]
	J = np.vstack((JR,Jp,Jz)).T

	# thin shade
	errfac1 = 0.1
	errfac2 = 0.01
	sprfac = 5

	if(with_boxes):
		thinJRerr=1./disc_dff(J,thin_sigmaR,kappa)
		thinJzerr=1./disc_dff(J,thin_sigmaz,nu)

		maxthinJ = np.sqrt(thinJRerr**2+thinJzerr**2)/Jp[0]

		a[0].fill_between(JRJz,errfac1*thinJRerr,errfac2*thinJRerr,where=JRJz<sprfac*maxthinJ,alpha=0.2,color=sns.color_palette()[0])
		a[1].fill_between(JRJz,errfac1*thinJzerr,errfac2*thinJzerr,where=JRJz<sprfac*maxthinJ,alpha=0.2,color=sns.color_palette()[0])
		a[0].annotate('Thin', xy=(0.0012,errfac2*thinJRerr[0]), fontsize=12,
		                horizontalalignment='left', verticalalignment='bottom',color=sns.color_palette()[0])
		a[1].annotate('Thin', xy=(0.0012,errfac2*thinJzerr[0]), fontsize=12,
		                horizontalalignment='left', verticalalignment='bottom',color=sns.color_palette()[0])

		thickJRerr=1./disc_dff(J,thick_sigmaR,kappa)
		thickJzerr=1./disc_dff(J,thick_sigmaz,nu)

		maxthickJ = np.sqrt(thickJRerr**2+thickJzerr**2)/Jp[0]

		a[0].fill_between(JRJz,errfac1*thickJRerr,errfac2*thickJRerr,where=JRJz<sprfac*maxthickJ,alpha=0.2,color=sns.color_palette()[1])
		a[1].fill_between(JRJz,errfac1*thickJzerr,errfac2*thickJzerr,where=JRJz<sprfac*maxthickJ,alpha=0.2,color=sns.color_palette()[1])

		a[0].annotate('Thick', xy=(0.0012,errfac1*thickJRerr[0]*0.9), fontsize=12,
		                horizontalalignment='left', verticalalignment='top',color=sns.color_palette()[1])
		a[1].annotate('Thick', xy=(0.0012,errfac1*thickJzerr[0]*0.9), fontsize=12,
		                horizontalalignment='left', verticalalignment='top',color=sns.color_palette()[1])

		haloJRerr=1./halo_dff(J,alpha0,alphainf,f,delta,v0,r0,rad=True)
		haloJzerr=1./halo_dff(J,alpha0,alphainf,f,delta,v0,r0,rad=False)

		a[0].fill_between(JRJz,errfac1*haloJRerr,errfac2*haloJRerr,alpha=0.4,color=sns.color_palette()[4])
		a[1].fill_between(JRJz,errfac1*haloJzerr,errfac2*haloJzerr,alpha=0.4,color=sns.color_palette()[4])
		a[0].annotate('Halo', xy=(0.0012,errfac2*haloJRerr[0]), fontsize=12,
		                horizontalalignment='left', verticalalignment='bottom',color=sns.color_palette()[4])
		a[1].annotate('Halo', xy=(0.0012,errfac2*haloJzerr[0]), fontsize=12,
		                horizontalalignment='left', verticalalignment='bottom',color=sns.color_palette()[4])

		streamJRerr=(results.T[5]-results.T[4])/np.pi
		streamLzerr = results.T[4]/np.pi
		streamJzerr=2.*results.T[6]/np.pi
		sigmav1=5.
		sigmav2=0.5

		a[0].set_xlim(0.001,0.6) #JRJz[-1])
		a[1].set_xlim(0.001,0.6) #JRJz[-1])
		a[0].fill_between(JRJz,sigmav1*streamJRerr,sigmav2*streamJRerr,where=JRJz>0.01,alpha=0.2,color=sns.color_palette()[2])
		a[1].fill_between(JRJz,sigmav1*streamJzerr,sigmav2*streamJzerr,where=JRJz>0.01,alpha=0.2,color=sns.color_palette()[2])
		if(JRJz[-1]>0.01):
			a[0].annotate('Stream', xy=(JRJz[JRJz>0.01][0],sigmav1*1.3*streamJRerr[JRJz>0.01][0]), fontsize=12,
		                horizontalalignment='left', verticalalignment='top',color=sns.color_palette()[2],rotation=20)
			a[1].annotate('Stream', xy=(JRJz[JRJz>0.01][0],sigmav1*1.3*streamJzerr[JRJz>0.01][0]), fontsize=12,
		                horizontalalignment='left', verticalalignment='top',color=sns.color_palette()[2],rotation=20)

	offsets=[1.,1.5,1.,1.15,0.85,1.,1.,1.]
	offsets2=[1.,1.,1.,1.1,0.9,1.,1.,1.]

	offsetsR=[1.,1.,1.4,0.75,1.2,1.,1.,1.]
	offsetsz=[1.,0.97,1.5,0.8,0.8,1.,1.,1.]

	for N,i in enumerate([5,6,0,4,7,1,3,2]):
		l1, = a[0].plot(JRJz,results.T[i*8+11],label=labels[i],lw=.5,ls=dashed[i],color=colors[i])
		a[0].annotate(labels[i],xy=(MaxX*1.012,results.T[i*8+11][-1]*offsetsR[N]),horizontalalignment='left', verticalalignment=align[i],color=colors[i],fontsize=5)
		l2, = a[1].plot(JRJz,results.T[i*8+12],lw=0.5,ls=dashed[i],color=colors[i])
		a[1].annotate(labels[i],xy=(MaxX*1.012,results.T[i*8+12][-1]*offsetsz[N]), horizontalalignment='left', verticalalignment=align[i],color=colors[i],fontsize=5)
		if(i==1 or i==2 or i==3):
			l1.set_dashes((2,1))
			l2.set_dashes((2,1))

	plt.setp(a[0].get_xticklabels(),visible=False)
	plt.setp(a[0].get_yticklabels()[1],visible=False)

	if(with_res):
		draw_res_lines(a,res_line,1.2e-4)

	a[1].set_xlabel(r'$(J_R+J_z)/|J_\phi|$')
	a[0].set_ylabel(r'$\sqrt{\langle(\Delta J_R)^2\rangle}/\mathrm{kpc\,km\,s}^{-1}$')
	a[1].set_ylabel(r'$\sqrt{\langle(\Delta J_z)^2\rangle}/\mathrm{kpc\,km\,s}^{-1}$')
	a[0].semilogx()
	a[1].semilogx()
	a[0].semilogy()
	a[1].semilogy()
	leg=a[0].legend(handlelength=0.8, scatterpoints=1, numpoints=1,frameon=False,ncol=4,loc='lower center', bbox_to_anchor=(0.5, 1.0))
	# set the linewidth of each legend object
	for legobj in leg.legendHandles:
	    legobj.set_linewidth(2.0)
	a[0].set_ylim(0.0001,1000.)
	a[1].set_ylim(0.0001,1000.)
	plt.savefig(results_file+'.acc.pdf',bbox_inches='tight')

	fig,a = plt.subplots(1,1,figsize=(4.,4))

	offsets=np.ones(8)
	offsets[2]=1.2
	offsets[3]=0.8
	offsets[1]=1.3
	offsets[-1]=1.1

	for N,i in enumerate([5,6,0,4,7,1,3,2]):
		n=i
		if(i>5):
			n=i-6
		l1, = a.plot(JRJz,results.T[75+i]/10**9,label=labels[i],lw=.5,ls=dashed[i],color=colors[i])
		if(i==1 or i==2 or i==3):
			l1.set_dashes((2,1))
			l2.set_dashes((2,1))
		plt.annotate(labels[i],xy=(MaxX*1.012,offsets[N]*results.T[75+i][-1]/10**9),horizontalalignment='left', verticalalignment=align[i],color=colors[i],fontsize=5)
	a.set_xlim(0.001,0.6) #JRJz[-1])
	a.set_ylim(1e-6,1.)
	a.set_xlabel(r'$(J_R+J_z)/|J_\phi|$')
	a.set_ylabel(r'$t/\mathrm{s}$')
	a.semilogx()
	a.semilogy()
	leg=a.legend(handlelength=0.8, scatterpoints=1, numpoints=1,frameon=False,ncol=4,loc='lower center', bbox_to_anchor=(0.5, 1.0))
	# set the linewidth of each legend object
	for legobj in leg.legendHandles:
	    legobj.set_linewidth(2.0)
	plt.savefig(results_file+'.times.pdf',bbox_inches='tight')

	fig,a = plt.subplots(3,1,figsize=(4.,8.))

	plt.subplots_adjust(hspace=0.)
	colors = [sns.color_palette()[2],sns.color_palette()[0],sns.color_palette()[1],sns.color_palette("Set2", 10)[3],sns.color_palette()[3],sns.color_palette()[5],sns.color_palette()[4],sns.color_palette("Set1",5,desat=.6)[4]]


	dashed=['-','--','--','--','-','-','-','-']

	JRJz=results.T[3]
	JR = results.T[0]
	Jz = results.T[2]
	Jp = results.T[1]
	J = np.vstack((JR,Jp,Jz)).T

	a[0].set_xlim(0.001,0.6) #JRJz[-1])
	a[1].set_xlim(0.001,0.6) #JRJz[-1])
	a[2].set_xlim(0.001,0.6) #JRJz[-1])

	offsetsR=[1.,1.,1.5,0.9,0.8,1.,1.,1.]
	offsetsp=[0.85,1.,1.7,0.9,0.9,1.,1.,1.]
	offsetsz=[1.4,1.,1.,1.,1.,1.,1.,1.]

	a[0].plot(JRJz,kms2kpcGyr*results.T[7],c='k')
	a[1].plot(JRJz,kms2kpcGyr*results.T[8],c='k')
	a[2].plot(JRJz,kms2kpcGyr*results.T[9],c='k')

	a[0].annotate(r'$\Omega_R$',xy=(0.2,kms2kpcGyr*results.T[7][370]*0.8),horizontalalignment='right', verticalalignment='top',fontsize=12)
	a[1].annotate(r'$\Omega_\phi$',xy=(0.2,kms2kpcGyr*results.T[8][370]*0.8),horizontalalignment='right', verticalalignment='top',fontsize=12)
	a[2].annotate(r'$\Omega_z$',xy=(0.2,kms2kpcGyr*results.T[9][370]*0.7),horizontalalignment='right', verticalalignment='top',fontsize=12)

	for N,i in enumerate([5,6,0,4,7,1,3,2]):
		if(i==3):
			continue
		l1, = a[0].plot(JRJz,kms2kpcGyr*results.T[i*8+16],label=labels[i],lw=.5,ls=dashed[i],color=colors[i])
		a[0].annotate(labels[i],xy=(MaxX*1.012,kms2kpcGyr*results.T[i*8+16][-1]*offsetsR[N]),horizontalalignment='left', verticalalignment=align[i],color=colors[i],fontsize=5)
		l2, = a[1].plot(JRJz,kms2kpcGyr*results.T[i*8+17],lw=0.5,ls=dashed[i],color=colors[i])
		a[1].annotate(labels[i],xy=(MaxX*1.012,kms2kpcGyr*results.T[i*8+17][-1]*offsetsp[N]), horizontalalignment='left', verticalalignment=align[i],color=colors[i],fontsize=5)
		l3, = a[2].plot(JRJz,kms2kpcGyr*results.T[i*8+18],label=labels[i],lw=.5,ls=dashed[i],color=colors[i])
		a[2].annotate(labels[i],xy=(MaxX*1.012,kms2kpcGyr*results.T[i*8+18][-1]*offsetsz[N]),horizontalalignment='left', verticalalignment=align[i],color=colors[i],fontsize=5)
		if(i==1 or i==2 or i==3):
			l1.set_dashes((2,1))
			l2.set_dashes((2,1))
			l3.set_dashes((2,1))

	plt.setp(a[0].get_xticklabels(),visible=False)
	plt.setp(a[1].get_xticklabels(),visible=False)
	plt.setp(a[0].get_yticklabels()[1],visible=False)
	plt.setp(a[1].get_yticklabels()[1],visible=False)

	a[2].set_xlabel(r'$(J_R+J_z)/|J_\phi|$')
	a[0].set_ylabel(r'$\sqrt{\langle(\Delta \Omega_R)^2\rangle}/\mathrm{Gyr}^{-1}$')
	a[1].set_ylabel(r'$\sqrt{\langle(\Delta \Omega_\phi)^2\rangle}/\mathrm{Gyr}^{-1}$')
	a[2].set_ylabel(r'$\sqrt{\langle(\Delta \Omega_z)^2\rangle}/\mathrm{Gyr}^{-1}$')

	if(with_res):
		draw_res_lines(a,res_line,1.2e-6)
	if(with_boxes):
		D = np.reshape(results.T[-9:].T,(len(results),3,3))
		Dinv = np.linalg.inv(D)
		sJ = np.array([np.diag([i,j,k]) for i,j,k in zip(streamJRerr,streamLzerr,streamJzerr)])
		Omg = np.array([np.dot(i,np.dot(j,i)) for i,j in zip(Dinv,sJ)])
		Omg = np.linalg.inv(Omg)

		streamOmRerr=np.sqrt(Omg[:,0,0])
		streamOmperr=np.sqrt(Omg[:,1,1])
		streamOmzerr=np.sqrt(Omg[:,2,2])

		a[0].fill_between(JRJz,sigmav1*streamOmRerr,sigmav2*streamOmRerr,where=(JRJz>0.01)*(JRJz<0.38),alpha=0.2,color=sns.color_palette()[2])
		a[1].fill_between(JRJz,sigmav1*streamOmperr,sigmav2*streamOmperr,where=(JRJz>0.01)*(JRJz<0.38),alpha=0.2,color=sns.color_palette()[2])
		a[2].fill_between(JRJz,sigmav1*streamOmzerr,sigmav2*streamOmzerr,where=(JRJz>0.01)*(JRJz<0.38),alpha=0.2,color=sns.color_palette()[2])
	# if(JRJz[-1]>0.01):
	# 	a[0].annotate('Stream', xy=(JRJz[1:][JRJz[1:]>0.01][0],sigmav1*2.7*streamOmRerr[JRJz[1:]>0.01][0]), fontsize=12,
	#                 horizontalalignment='left', verticalalignment='top',color=sns.color_palette()[4],rotation=10)
	# 	a[1].annotate('Stream', xy=(JRJz[1:][JRJz[1:]>0.01][0]*1.05,sigmav1*3.5*streamOmperr[JRJz[1:]>0.01][0]), fontsize=12,
	#                 horizontalalignment='left', verticalalignment='top',color=sns.color_palette()[4],rotation=4)
	# 	a[2].annotate('Stream', xy=(JRJz[1:][JRJz[1:]>0.01][0],sigmav1*2.6*streamOmzerr[JRJz[1:]>0.01][0]), fontsize=12,
	#                 horizontalalignment='left', verticalalignment='top',color=sns.color_palette()[4],rotation=-10)

	a[0].semilogx()
	a[1].semilogx()
	a[2].semilogx()
	a[0].semilogy()
	a[1].semilogy()
	a[2].semilogy()
	leg=a[0].legend(handlelength=0.8, scatterpoints=1, numpoints=1,frameon=False,ncol=4,loc='lower center', bbox_to_anchor=(0.5, 1.0))
	# set the linewidth of each legend object
	for legobj in leg.legendHandles:
	    legobj.set_linewidth(2.0)
	# a[0].set_ylim(0.0001,1000.)
	# a[1].set_ylim(0.0001,1000.)
	# a[2].set_ylim(0.0001,1000.)


	plt.savefig(results_file+'.freq.pdf',bbox_inches='tight')

	fig,a = plt.subplots(3,1,figsize=(4.,8.))

	plt.subplots_adjust(hspace=0.)
	colors = [sns.color_palette()[2],sns.color_palette()[0],sns.color_palette()[1],sns.color_palette("Set2", 10)[3],sns.color_palette()[3],sns.color_palette()[5],sns.color_palette()[4],sns.color_palette("Set1",5,desat=.6)[4]]


	dashed=['-','--','--','--','-','-','-','-']

	JRJz=results.T[3]
	JR = results.T[0]
	Jz = results.T[2]
	Jp = results.T[1]
	J = np.vstack((JR,Jp,Jz)).T

	a[0].set_xlim(0.001,0.6)
	a[1].set_xlim(0.001,0.6)
	a[2].set_xlim(0.001,0.6)

	offsets=[1.,0.85,1.3,0.85,1.3,1.,1.,1.]
	offsets2=[1.,1.,1.1,0.85,0.9,1.,1.,1.]
	offsets3=[1.,1.,1.2,0.85,0.9,1.,1.,1.]

	for N,i in enumerate([5,6,0,4,7,1,3,2]):
		if(i==3):
			continue
		l1, = a[0].plot(JRJz,results.T[i*8+13]/np.pi,label=labels[i],lw=.5,ls=dashed[i],color=colors[i])
		a[0].annotate(labels[i],xy=(MaxX*1.012,results.T[i*8+13][-1]*offsets2[N]/np.pi),horizontalalignment='left', verticalalignment=align[i],color=colors[i],fontsize=5)
		l2, = a[1].plot(JRJz,results.T[i*8+14]/np.pi,lw=0.5,ls=dashed[i],color=colors[i])
		a[1].annotate(labels[i],xy=(MaxX*1.012,results.T[i*8+14][-1]*offsets[N]/np.pi), horizontalalignment='left', verticalalignment=align[i],color=colors[i],fontsize=5)
		l3, = a[2].plot(JRJz,results.T[i*8+15]/np.pi,label=labels[i],lw=.5,ls=dashed[i],color=colors[i])
		a[2].annotate(labels[i],xy=(MaxX*1.012,results.T[i*8+15][-1]*offsets3[N]/np.pi),horizontalalignment='left', verticalalignment=align[i],color=colors[i],fontsize=5)
		if(i==1 or i==2 or i==3):
			l1.set_dashes((2,1))
			l2.set_dashes((2,1))
			l3.set_dashes((2,1))
	if(with_boxes):
		streamThRerr=1./(results.T[5]-results.T[4])
		streamThperr=1./results.T[4]
		streamThzerr=.5/results.T[6]
		r1=0.1
		r2=0.01

		a[0].fill_between(JRJz,r1*streamThRerr,r2*streamThRerr,where=JRJz>0.01,alpha=0.2,color=sns.color_palette()[2])
		a[1].fill_between(JRJz,r1*streamThperr,r2*streamThperr,where=JRJz>0.01,alpha=0.2,color=sns.color_palette()[2])
		a[2].fill_between(JRJz,r1*streamThzerr,r2*streamThzerr,where=JRJz>0.01,alpha=0.2,color=sns.color_palette()[2])
	# if(JRJz[-1]>0.01):
	# 	a[0].annotate('Stream', xy=(JRJz[JRJz>0.01][0],r1*0.9*streamThRerr[JRJz>0.01][0]), fontsize=12,
	#                 horizontalalignment='left', verticalalignment='top',color=sns.color_palette()[4],rotation=-15)
	# 	a[1].annotate('Stream', xy=(JRJz[JRJz>0.01][0]*1.05,r1*1.*streamThperr[JRJz>0.01][0]), fontsize=12,
	#                 horizontalalignment='left', verticalalignment='top',color=sns.color_palette()[4],rotation=0)
	# 	a[2].annotate('Stream', xy=(JRJz[JRJz>0.01][0],r1*0.45*streamThzerr[JRJz>0.01][0]), fontsize=12,
	#                 horizontalalignment='left', verticalalignment='top',color=sns.color_palette()[4],rotation=-20)


	plt.setp(a[0].get_xticklabels(),visible=False)
	plt.setp(a[1].get_xticklabels(),visible=False)
	plt.setp(a[0].get_yticklabels()[1],visible=False)
	plt.setp(a[1].get_yticklabels()[1],visible=False)

	a[2].set_xlabel(r'$(J_R+J_z)/|J_\phi|$')
	a[0].set_ylabel(r'$\sqrt{\langle(\Delta \theta_R)^2\rangle}/\pi$')
	a[1].set_ylabel(r'$\sqrt{\langle(\Delta \theta_\phi)^2\rangle}/\pi$')
	a[2].set_ylabel(r'$\sqrt{\langle(\Delta \theta_z)^2\rangle}/\pi$')

	if(with_res):
		draw_res_lines(a,res_line,0.3)

	a[0].semilogx()
	a[1].semilogx()
	a[2].semilogx()
	a[0].semilogy()
	a[1].semilogy()
	a[2].semilogy()
	leg=a[0].legend(handlelength=0.8, scatterpoints=1, numpoints=1,frameon=False,ncol=4,loc='lower center', bbox_to_anchor=(0.5, 1.0))
	# set the linewidth of each legend object
	for legobj in leg.legendHandles:
	    legobj.set_linewidth(2.0)
	# a[0].set_ylim(0.0001,1000.)
	# a[1].set_ylim(1e-6,.1)
	# a[2].set_ylim(0.0001,1000.)


	plt.savefig(results_file+'.ang.pdf',bbox_inches='tight')

def genfunc_converg(data,labels):
	plt.clf()
	f = plt.figure(figsize=(3.32,2.))
	markers = ['o','h','s','^']
	Text = [r'$2$',
			r'$4$',
			r'$6$',
			r'$8$',
			r'$10$',
			r'$12$',
			r'$14$',
			r'$16$',
			r'$18$']
	for n,(d,l) in enumerate(zip(data,labels)):
		dat = np.genfromtxt(d)
		plt.plot(dat.T[3]/1.e9,dat.T[1],color=sns.color_palette()[n])
		plt.plot(dat.T[3]/1.e9,dat.T[2],color=sns.color_palette()[n])
		plt.scatter(dat.T[3]/1.e9,dat.T[1],
		            color=sns.color_palette()[n],marker=markers[n],
		            label=labels[n])
		plt.scatter(dat.T[3]/1.e9,dat.T[2],
		            color=sns.color_palette()[n],marker=markers[n])
		if(n==2):
			for i, txt in enumerate(Text):
			    if(i==0):
			    	plt.annotate(txt, (dat.T[3][i]/1.e9*1.07,dat.T[2][i]*0.7),fontsize=6)
			    else:
			    	plt.annotate(txt, (dat.T[3][i]/1.e9*1.07,dat.T[2][i]),fontsize=6)
	plt.annotate(r'$N_\mathrm{samp}=2400,\,N_T=24$',xy=(0.95,0.95), xycoords='axes fraction', horizontalalignment='right', verticalalignment='top')
	plt.xlabel(r'$t/s$')
	plt.ylabel(r'$\Delta J_i/\mathrm{kpc\,km\,s}^{-1}$')
	plt.semilogy()
	plt.semilogx()
	plt.legend(handlelength=1, scatterpoints=1, numpoints=1,frameon=False,ncol=4,loc='lower center', bbox_to_anchor=(0.5, 1.),fontsize=10)
	plt.savefig('genfunc_converg.pdf',bbox_inches='tight')

import tabulate

def make_table(files,names,dontplot,Jranges):
	actmethods=['PAA','SAA','Fudge v1','Fudge v2','Fit','ItTorus','AvGenfunc','Genfunc']
	acts=[[0 for y in range(9)]]
	acts[0][0]='Actions'
	spreadsnon=[[0 for x in range(9)] for y in range(5)]
	spreadscon=[[0 for x in range(9)] for y in range(3)]
	for i in range(len(actmethods)):
		if(i>4):
			spreadscon[i-5][0]=actmethods[i]
		else:
			spreadsnon[i][0]=actmethods[i]
	for N,(f,n,d) in enumerate(zip(files,names,dontplot)):
		r = method_comparison(f,n,dontplot=d,rangess=5.,Jranges=Jranges[N])
		acts[0][1+2*N]=r[0][0]
		acts[0][2+2*N]=r[0][1]
		for i in range(len(r)-1):
			if(i>4):
				spreadscon[i-5][1+2*N]=r[i+1][0]
				spreadscon[i-5][2+2*N]=r[i+1][1]
			else:
				spreadsnon[i][1+2*N]=r[i+1][0]
				spreadsnon[i][2+2*N]=r[i+1][1]
	print acts
	actstab= tabulate.tabulate(acts,names,tablefmt='latex',floatfmt='.2f')
	tablenon= tabulate.tabulate(spreadsnon,tablefmt='latex',floatfmt='.1g')
	tablecon= tabulate.tabulate(spreadscon,tablefmt='latex',floatfmt='.1g')
	outfile=open('acttable','w')
	outfile.write(actstab)
	outfile.write(tablenon)
	outfile.write(tablecon)
	outfile.close

if __name__=="__main__":
    time_plot('many_tori_output.dat', with_res=True, with_boxes=True)
    time_plot('many_tori_output_james.dat',with_res=True,with_boxes=True)
    make_table(['thin','thick','halo','stream'],['Thin','Thick','Halo','Stream'],[[],[],[5],[5]],[None,None,None,[[290.,340.],[540.,580.]]])
    genfunc_converg(['thin_gc','thick_gc','halo_gc','stream_gc'],['Thin','Thick','Halo','Stream'])
