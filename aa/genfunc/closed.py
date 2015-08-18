import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.axes3d import Axes3D

col0 = "black"
colA = "#377eb8"
colB = "#4daf4a"
colC = "#e41a1c"

fig_width_pt = 480.0            # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0 / 72.27       # Convert pt to inch
golden_mean = (np.sqrt(5) - 1.0) / 2.0         # Aesthetic ratio
fig_width = fig_width_pt * inches_per_pt  # width in inches
fig_height = fig_width * golden_mean      # height in inches
fig_size = [fig_width,0.7 * fig_height]

params = {'backend': 'Cairo','axes.labelsize': 10,'text.fontsize': 10,'legend.fontsize': 10,
          'xtick.labelsize': 10,'ytick.labelsize': 10,'text.usetex': False,'figure.figsize': fig_size,
          'font.family':'Times New Roman','font.serif':'Times'}
rcParams.update(params)


def plot_closed_short():
    plt.figure(figsize=[3.32,3.32])

    Delta = np.genfromtxt("Delta1.dat")
    plt.plot(Delta.T[0],np.log10(Delta.T[1]),'k')
    plt.xlabel(r'$\Delta_1/{\rm kpc}$')
    plt.ylabel(r'$\log10(\Delta J_\mu/{\rm kpc\,km\,s}^{-1})$')
    ax=plt.gca()
    plt.text(0, 1.02, 'Closed short-axis loop', horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)

    a2 = plt.axes([0.55,0.25,0.3,0.3])
    Orb = np.genfromtxt("orbit_Delta1.dat")
    a2.plot(Orb.T[0],Orb.T[1],'k',zorder=0)
    a2.arrow(0.,0.,0.,Orb[1][1] - 0.5,color='r',linewidth=0.2,head_width=0.5, head_length=0.5,zorder=1)
    a2.arrow(0.,0.,Orb[1][1] - 0.5,0.,color='r',linewidth=0.2,head_width=0.5, head_length=0.5,zorder=1)
    a2.set_aspect('equal')
    a2.set_xlabel(r'$x/{\rm kpc}$')
    a2.set_ylabel(r'$y/{\rm kpc}$')

    plt.savefig('Delta1.pdf',bbox_inches='tight')
    plt.clf()


def plot_closed_long():
    plt.figure(figsize=[3.32,3.32])

    Delta = np.genfromtxt("Delta2.dat")
    Delta1 = np.genfromtxt("Delta2_loweralpha.dat")
    plt.plot(Delta.T[0],np.log10(Delta.T[1]),'k',label=r'$\alpha=-80\,{\rm kpc}^2$')
    plt.plot(Delta1.T[0],np.log10(Delta1.T[1]),'rx',label=r'$\alpha=-20\,{\rm kpc}^2$')
    ax=plt.gca()
    plt.text(0, 1.02, 'Closed long-axis loop', horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)

    lg = plt.legend(loc=2)
    lg._drawFrame=False
    plt.xlabel(r'$\Delta_2/{\rm kpc}$')
    plt.ylabel(r'$\log_{10}(\Delta J_\nu/{\rm kpc\,km\,s}^{-1})$')

    a2 = plt.axes([0.6,0.25,0.3,0.3])
    Orb = np.genfromtxt("orbit_Delta2.dat")
    a2.plot(Orb.T[1],Orb.T[2],'k',zorder=0)
    a2.arrow(0.,0.,0.,Orb[1][1] - 0.5,color='r',linewidth=0.2,head_width=0.5, head_length=0.5,zorder=1)
    a2.arrow(0.,0.,Orb[1][1] - 0.5,0.,color='r',linewidth=0.2,head_width=0.5, head_length=0.5,zorder=1)
    a2.set_aspect('equal')
    a2.set_xlabel(r'$y/{\rm kpc}$')
    a2.set_ylabel(r'$z/{\rm kpc}$')

    plt.savefig('Delta2.pdf',bbox_inches='tight')
    plt.clf()


def Delta_withEnergy():
    plt.figure(figsize=[3.32,3.32 * 0.7])
    Delta = np.genfromtxt("Delta_withE.dat")
    plt.plot(Delta.T[0] / 10000.,np.sqrt(Delta.T[1] - Delta.T[2]),'k',label=r'$\Delta_1$')
    plt.plot(Delta.T[0] / 10000.,np.sqrt(-1. - Delta.T[1]),color=colC,linestyle='--',label=r'$\Delta_2$')
    plt.axvline(-8.40295,color=colA,linestyle=':')
    lg = plt.legend(loc=2,handlelength=2.5)
    lg._drawFrame=False
    plt.xlabel(r'$E/{\rm (100\,km\,s}^{-1})^2$')
    plt.ylabel(r'$\Delta_i/{\rm kpc}$')
    plt.savefig('Delta_withE.pdf',bbox_inches='tight')
    plt.clf()


def DeltaCheck():

    jet = plt.get_cmap('jet')
    # c = colors.ColorConverter().to_rgb
    # rvb = make_colormap([c("#91bfdb"), c("#ffffbf"), 1. / 2., c("#ffffbf"),c("#fc8d59")])
    cNorm = colors.Normalize(vmin=0.,vmax=np.pi / 2.)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    scalarMap._A=[]
    f,a = plt.subplots(2,1,figsize=[3.32,3.32 * 0.7])
    plt.subplots_adjust(hspace=0.)
    Delta = np.genfromtxt("delta_fiddle_lower.dat")
    a[0].scatter(Delta.T[0],np.sqrt(Delta.T[3] - Delta.T[2]),c=Delta.T[1],s=4.3,edgecolor="None",zorder=1)
    a[0].axhline(np.sqrt(-13.3317 + 20.1094),color='k',linewidth=1,zorder=0)
    a[1].scatter(Delta.T[0],np.sqrt(-1. - Delta.T[3]),c=Delta.T[1],label=r'$\Delta_2/{\rm kpc}$',s=4.3,edgecolor="None",zorder=1)
    a[1].axhline(np.sqrt(-1 + 13.3317),linewidth=1,color='k',zorder=1)
    a[0].set_xticks([])
    a[1].set_yticks([0.,1.,2.,3.,4.,5.])
    a[0].set_xlim(0.,10.)
    a[0].set_ylim(0.,5.)
    a[1].set_xlim(0.,10.)
    a[1].set_ylim(0.,6.)
    plt.xlabel(r'$y/{\rm kpc}$')
    a[0].set_ylabel(r'$\Delta_1/{\rm kpc}$')
    a[1].set_ylabel(r'$\Delta_2/{\rm kpc}$')
    cbar_ax = f.add_axes([0.95, 0.15, 0.05, 0.7])
    ccc = f.colorbar(scalarMap,cax=cbar_ax)
    ccc.set_label(r'$\theta/{\rm rad}$')
    plt.savefig('DeltaCheck.pdf',bbox_inches='tight')
    plt.clf()


def orbit_accuracy_log10fractional():
    jet = plt.get_cmap('jet')
    # c = colors.ColorConverter().to_rgb
    # rvb = make_colormap([c("#91bfdb"), c("#ffffbf"), 1. / 2., c("#ffffbf"),c("#fc8d59")])
    cNorm = colors.Normalize(vmin=0., vmax=np.pi / 2.)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    scalarMap._A=[]
    f,a = plt.subplots(3,1,figsize=[3.32,3.32 * 1.4])
    plt.subplots_adjust(hspace=0.)
    Delta = np.genfromtxt("orbits.dat")
    a[0].scatter(Delta.T[0],np.log10(Delta.T[2] / Delta.T[5]),c=Delta.T[1],s=4.3,edgecolor="None",zorder=1)
    a[1].scatter(Delta.T[0],np.log10(Delta.T[3] / Delta.T[6]),c=Delta.T[1],s=4.3,edgecolor="None",zorder=1)
    a[2].scatter(Delta.T[0],np.log10(Delta.T[4] / Delta.T[7]),c=Delta.T[1],s=4.3,edgecolor="None",zorder=1)
    a[0].set_xticks([])
    yticks = a[1].yaxis.get_major_ticks()
    yticks[-1].label1.set_visible(False)
    a[1].set_xticks([])
    yticks = a[2].yaxis.get_major_ticks()
    yticks[-1].label1.set_visible(False)
    plt.xlabel(r'$y/{\rm kpc}$')
    a[0].set_ylabel(r'$\log_{10}(\Delta J_\lambda/J_\lambda)$')
    a[1].set_ylabel(r'$\log_{10}(\Delta J_\mu/J_\mu)$')
    a[2].set_ylabel(r'$\log_{10}(\Delta J_\nu/J_\nu)$')
    cbar_ax = f.add_axes([0.15, 0.97, 0.7, 0.03])
    ccc = f.colorbar(scalarMap,cax=cbar_ax,orientation='horizontal')
    ccc.set_label(r'$\theta/{\rm rad}$')
    ccc.ax.xaxis.set_ticks_position('top')
    ccc.ax.xaxis.set_label_position('top')
    plt.savefig('orbit_accuracy_log10.pdf',bbox_inches='tight')
    plt.clf()


def orbit_accuracy_absolute():
    f,a = plt.subplots(3,1,figsize=[3.32,3.32 * 1.65907])
    plt.subplots_adjust(hspace=0.)
    Delta = np.genfromtxt("orbits.dat")
    a[0].scatter(Delta.T[0],np.log10(Delta.T[5]),c=colC,s=4.3,edgecolor="None",zorder=0,label=r"$J_i$")
    a[1].scatter(Delta.T[0],np.log10(Delta.T[6]),c=colC,s=4.3,edgecolor="None",zorder=0)
    a[2].scatter(Delta.T[0],np.log10(Delta.T[7]),c=colC,s=4.3,edgecolor="None",zorder=0)
    a[0].scatter(Delta.T[0],np.log10(Delta.T[2]),c='k',marker='^',s=4.3,edgecolor="None",zorder=1,label=r"$\Delta J_i$")
    a[1].scatter(Delta.T[0],np.log10(Delta.T[3]),c='k',marker='^',s=4.3,edgecolor="None",zorder=1)
    a[2].scatter(Delta.T[0],np.log10(Delta.T[4]),c='k',marker='^',s=4.3,edgecolor="None",zorder=1)

    a[0].set_xticks([])
    yticks = a[1].yaxis.get_major_ticks()
    yticks[-1].label1.set_visible(False)
    a[1].set_xticks([])
    yticks = a[2].yaxis.get_major_ticks()
    yticks[-1].label1.set_visible(False)
    plt.xlabel(r'$y/{\rm kpc}$')
    a[0].set_ylabel(r'$\log_{10}(J_\lambda/{\rm kpc\,km\,s}^{-1})$')
    a[1].set_ylabel(r'$\log_{10}(J_\mu/{\rm kpc\,km\,s}^{-1})$')
    a[2].set_ylabel(r'$\log_{10}(J_\nu/{\rm kpc\,km\,s}^{-1})$')
    lg = a[0].legend(loc='lower center', bbox_to_anchor=(0.5, 1.05),ncol=2)
    lg._drawFrame=False
    plt.savefig('orbit_accuracy.pdf',bbox_inches='tight')
    plt.clf()


def orbit_accuracy_sum():
    f,a = plt.subplots(3,1,figsize=[3.32,3.32 * 1.4])
    plt.subplots_adjust(hspace=0.)
    Delta = np.genfromtxt("orbits.dat")
    a[0].scatter(Delta.T[5],np.log10(Delta.T[2]),c='k',marker='^',s=4.3,edgecolor="None",zorder=1)
    a[1].scatter(Delta.T[5],np.log10(Delta.T[3]),c='k',marker='^',s=4.3,edgecolor="None",zorder=1)
    a[2].scatter(Delta.T[5],np.log10(Delta.T[4]),c='k',marker='^',s=4.3,edgecolor="None",zorder=1)
    a[0].set_xticks([])
    yticks = a[1].yaxis.get_major_ticks()
    yticks[-1].label1.set_visible(False)
    a[1].set_xticks([])
    yticks = a[2].yaxis.get_major_ticks()
    yticks[-1].label1.set_visible(False)
    plt.xlabel(r'$(J_\lambda)/{\rm kpc\,km\,s}^{-1}$')
    a[0].set_ylabel(r'$\log_{10}(J_\lambda/{\rm kpc\,km\,s}^{-1})$')
    a[1].set_ylabel(r'$\log_{10}(J_\mu/{\rm kpc\,km\,s}^{-1})$')
    a[2].set_ylabel(r'$\log_{10}(J_\nu/{\rm kpc\,km\,s}^{-1})$')
    plt.savefig('orbit_accuracy_sum.pdf',bbox_inches='tight')
    plt.clf()


def testLoop(infile,truth,outFile,factor,MaxN,text2):

    f,a = plt.subplots(2,3,figsize=[6.64,3.32 * 1.4])
    plt.subplots_adjust(hspace=0.3,wspace=0.4)
    Delta = np.genfromtxt(infile)
    color1 = "#1f78b4"
    color2 = "#b2df8a"
    a[0,0].plot(Delta.T[0],Delta.T[1],c='k')
    a[0,1].plot(Delta.T[0],Delta.T[2],c='k')
    a[0,2].plot(Delta.T[1],Delta.T[2],c='k')
    a[0,0].text(0, 1.02, text2, horizontalalignment='left', verticalalignment='bottom', transform=a[0,0].transAxes, fontsize=10)
    for i in range(3):
        a[0,i].set_aspect(1.)

    a[0,0].set_xlim(-np.max(Delta.T[0]),np.max(Delta.T[0]))
    a[0,0].set_ylim(-np.max(Delta.T[1]),np.max(Delta.T[1]))
    a[0,1].set_xlim(-np.max(Delta.T[0]),np.max(Delta.T[0]))
    a[0,1].set_ylim(-np.max(Delta.T[2]),np.max(Delta.T[2]))
    a[0,2].set_xlim(-np.max(Delta.T[1]),np.max(Delta.T[1]))
    a[0,2].set_ylim(-np.max(Delta.T[2]),np.max(Delta.T[2]))

    a[0,0].set_xlabel(r'$x/{\rm kpc}$')
    a[0,0].set_ylabel(r'$y/{\rm kpc}$')
    a[0,1].set_xlabel(r'$x/{\rm kpc}$')
    a[0,1].set_ylabel(r'$z/{\rm kpc}$')
    a[0,2].set_xlabel(r'$y/{\rm kpc}$')
    a[0,2].set_ylabel(r'$z/{\rm kpc}$')

    if(MaxN):
        a[0,1].yaxis.set_major_locator(MaxNLocator(4))
        a[0,2].yaxis.set_major_locator(MaxNLocator(4))

    a[1,0].scatter(Delta.T[6] / factor,Delta.T[7],c=color1,s=2.3,edgecolor="None",zorder=1)
    a[1,0].scatter(Delta.T[15] / factor,Delta.T[16],c=color2,s=2.3,edgecolor="None",zorder=2)
    a[1,0].axhline(truth[1],color='k')
    a[1,0].axvline(truth[0],color='k')
    a[1,0].xaxis.set_major_locator(MaxNLocator(6))
    a[1,0].yaxis.set_major_locator(MaxNLocator(6))
    a[1,1].scatter(Delta.T[6] / factor,Delta.T[8],c=color1,s=2.3,edgecolor="None",zorder=1)
    a[1,1].scatter(Delta.T[15] / factor,Delta.T[17],c=color2,s=2.3,edgecolor="None",zorder=2)
    a[1,1].axhline(truth[2],color='k')
    a[1,1].axvline(truth[0],color='k')
    a[1,1].xaxis.set_major_locator(MaxNLocator(6))
    a[1,1].yaxis.set_major_locator(MaxNLocator(6))
    a[1,2].scatter(Delta.T[7],Delta.T[8],c=color1,s=2.3,edgecolor="None",zorder=1)
    a[1,2].scatter(Delta.T[16],Delta.T[17],c=color2,s=2.3,edgecolor="None",zorder=2)
    a[1,2].axhline(truth[2],color='k')
    a[1,2].axvline(truth[1],color='k')
    a[1,2].xaxis.set_major_locator(MaxNLocator(6))
    a[1,2].yaxis.set_major_locator(MaxNLocator(6))
    plt.xlabel(r'$y/{\rm kpc}$')
    a[1,0].set_xlabel(r'$J_\lambda/{\rm kpc\,km\,s}^{-1}$')
    a[1,0].set_ylabel(r'$J_\mu/{\rm kpc\,km\,s}^{-1}$')
    a[1,1].set_xlabel(r'$J_\lambda/{\rm kpc\,km\,s}^{-1}$')
    a[1,1].set_ylabel(r'$J_\nu/{\rm kpc\,km\,s}^{-1}$')
    a[1,2].set_xlabel(r'$J_\mu/{\rm kpc\,km\,s}^{-1}$')
    a[1,2].set_ylabel(r'$J_\nu/{\rm kpc\,km\,s}^{-1}$')
    plt.savefig(outFile,bbox_inches='tight',dpi=400)
    plt.clf()


def plot_SoS(infile1, infile2, outFile, axis, axis2, text):
    f = plt.figure(figsize=[2.32,2.32 * 0.7])
    c = colors.ColorConverter().to_rgb
    rvb = make_colormap([c("#91bfdb"), c("#ffffbf"), 1. / 2., c("#ffffbf"),c("#fc8d59")])
    trueSoS = np.genfromtxt(infile1)
    StackelSoS = np.genfromtxt(infile2)
    color_max = np.max(np.array([np.abs(StackelSoS[i * 201][0]) for i in range(0, len(StackelSoS) / 201, 1)]))
    color_min = np.min(np.array([np.abs(StackelSoS[i * 201][0]) for i in range(0, len(StackelSoS) / 201, 1)]))
    cNorm = colors.Normalize(vmin=color_min,vmax=color_max)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=rvb)
    scalarMap._A=[]
    StackelSoS = np.array([StackelSoS[i:i + 201] for i in range(0, len(StackelSoS), 201)])
    for i in StackelSoS:
        plt.plot(i[1:].T[0],i[1:].T[1],color=scalarMap.to_rgba(np.abs(i[0][0])),linewidth=0.7) #alpha = 0.3. linewidth = 0.1
    plt.plot(trueSoS.T[0],trueSoS.T[1],'k')
    plt.xlabel(axis + r'$/{\rm kpc}$')
    plt.ylabel(r'$p$' + axis2 + r'$/{\rm km\,s}^{-1}$')
    ax = plt.gca()
    plt.text(1, 1.02, text, horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes, fontsize=10)
    cbar_ax = f.add_axes([0.95, 0.15, 0.05, 0.7])
    ccc = f.colorbar(scalarMap,cax=cbar_ax)
    ccc.set_label(r'$|x_i|/{\rm kpc}$')
    plt.savefig(outFile,bbox_inches='tight',dpi=400)
    plt.clf()


def plot_jeans_solution(infile, outfile, index1, counter, text, text2):
    data = np.genfromtxt(infile)
    f,a = plt.subplots(2,1,figsize=[2.32,2.32])
    colors = ['k','r','b']
    NN=0
    for i in index1:
        a[0].plot(data.T[0],-data.T[i],color=colors[NN],linewidth=0,marker='o',markersize=2,markeredgecolor='none',zorder=0)
        a[0].plot(data.T[0],data.T[i + 3],color=colC,linestyle='-',linewidth=0.5,marker='o',markersize=0,markeredgecolor='none',zorder=1)
        a[1].plot(data.T[0],(data.T[i] + data.T[i + 3]) * 100. / data.T[i],color=colors[NN],linewidth=1,marker='o',markersize=2,markeredgecolor='none')
        NN+=1
    plt.subplots_adjust(hspace=0.)
    a[0].set_xticks([])
    yticks = a[1].yaxis.get_major_ticks()
    yticks[-1].label1.set_visible(False)
    a[0].set_ylabel(r'$-\rho\partial\Phi/\partial x_j$')  # /{\rm km^2s}^{-2}{\rm kpc}^{-4}$')
    a[1].set_xlabel(r'$x_i/{\rm kpc}$')
    a[1].set_ylabel('Percentage Error')
    a[0].set_yscale('log')
    a[0].set_xscale('log')
    a[1].set_xscale('log')
    a[0].set_xlim(0.2,40.)
    a[1].set_xlim(0.2,40.)
    plt.text(-0.1, 1.16, text2, horizontalalignment='left', verticalalignment='bottom', transform=a[0].transAxes, fontsize=10)
    plt.text(1, 1.02, text, horizontalalignment='right', verticalalignment='bottom', transform=a[0].transAxes, fontsize=10)
    plt.text(0.9, 0.75,r'$j=$' + str(counter[0]), horizontalalignment='right', verticalalignment='bottom', transform=a[0].transAxes, fontsize=10)
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()


def constituent_density(infile, outfile, text):
    data = np.genfromtxt(infile)
    plt.figure(figsize=[3.32,3.32 * 0.7])
    plt.plot(data.T[0],data.T[1] / data[0][4],color='b',linewidth=1,marker='o',markersize=2,markeredgecolor='none',label='Box')
    plt.plot(data.T[0],data.T[2] / data[0][4],color='r',linewidth=1,marker='o',markersize=2,markeredgecolor='none',label='Short-axis loop')
    plt.plot(data.T[0],data.T[3] / data[0][4],color='g',linewidth=1,marker='o',markersize=2,markeredgecolor='none',label='Long-axis loop')
    plt.plot(data.T[0],data.T[4] / data[0][4],color='k',linewidth=1,marker='o',markersize=2,markeredgecolor='none',label='Total')
    # plt.plot(data.T[0],0.5 / data.T[0],linewidth=1,linestyle=':',color='k')
    lg = plt.legend(loc=3)
    lg._drawFrame=False
    plt.ylabel(r'$\rho/\rho_{0.5{\rm kpc}}$')
    plt.xlabel(r'$x_i/{\rm kpc}$')
    plt.xscale('log')
    plt.yscale('log')
    ax=plt.gca()
    plt.xlim(0.2,40.)
    plt.ylim(1e-6,2.)
    plt.text(1, 1.02, text, horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes, fontsize=10)
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()


def twoDdensity(infile,infile2,text,text2,outfile):
    data = np.genfromtxt(infile)
    f,a = plt.subplots(1,3,figsize=[9.6,3.2])
    plt.subplots_adjust(wspace=0.3)
    X = np.reshape(data.T[0],(40,40))
    Y = np.reshape(data.T[1],(40,40))
    Z = np.reshape(np.log10(data.T[2]),(40,40))
    Z2 = np.reshape(np.log10(data.T[3]),(40,40))
    CC = a[0].contour(X,Y,Z,12)
    print CC.levels
    a[0].set_xlabel(r'$x/{\rm kpc}$')
    a[0].set_ylabel(r'$y/{\rm kpc}$')
    a[0].text(0, 1.02, text2, horizontalalignment='left', verticalalignment='bottom', transform=a[0].transAxes, fontsize=10)
    a[0].text(0.05, 0.85, r'$\log_{10}\rho(z=0)$', horizontalalignment='left', verticalalignment='bottom', transform=a[0].transAxes, fontsize=10)
    a[0].set_aspect('equal')
    CC = a[1].contour(X,Y,Z2,12)
    print CC.levels
    a[1].set_xlabel(r'$x/{\rm kpc}$')
    a[1].set_ylabel(r'$z/{\rm kpc}$')
    a[1].text(0.05, 0.85, r'$\log_{10}\rho(y=0)$', horizontalalignment='left', verticalalignment='bottom', transform=a[1].transAxes, fontsize=10)
    a[1].set_aspect('equal')

    data = np.genfromtxt(infile2)
    a[2].plot(data.T[0],data.T[1] / data[0][4],color=colA,linewidth=1,marker='o',markersize=2,markeredgecolor='none',label='Box')
    a[2].plot(data.T[0],data.T[2] / data[0][4],color=colC,linewidth=1,marker='o',markersize=2,markeredgecolor='none',label='Short-axis loop')
    a[2].plot(data.T[0],data.T[3] / data[0][4],color=colB,linewidth=1,marker='o',markersize=2,markeredgecolor='none',label='Long-axis loop')
    a[2].plot(data.T[0],data.T[4] / data[0][4],color='k',linewidth=1,marker='o',markersize=2,markeredgecolor='none',label='Total')
    # plt.plot(data.T[0],0.5 / data.T[0],linewidth=1,linestyle=':',color='k')
    lg = plt.legend(loc=3)
    lg._drawFrame=False
    a[2].set_ylabel(r'$\rho/\rho_{0.5{\rm kpc}}$')
    a[2].set_xlabel(r'$x_i/{\rm kpc}$')
    a[2].set_xscale('log')
    a[2].set_yscale('log')
    a[2].set_xlim(0.2,40.)
    a[2].set_ylim(1e-6,2.)
    a[2].text(1, 1.02, text, horizontalalignment='right', verticalalignment='bottom', transform=a[2].transAxes, fontsize=10)

    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()


def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return colors.LinearSegmentedColormap('CustomMap', cdict)

# from scipy.interpolate import griddata


def error_planes():
    f = plt.figure(figsize=[14,3.17])
    plt.subplots_adjust(wspace=0.,bottom=-0.1)
    Delta = np.genfromtxt("orbits_scaled.dat")
    az,el = 47, 25
    ax = f.add_subplot(1, 4, 1, projection='3d')
    c = colors.ColorConverter().to_rgb
    rvb = make_colormap([c(col0), c(colA), 1. / 3., c(colA), c(colB), 2. / 3., c(colB),c(colC)])
    cNorm = colors.Normalize(vmin=0.,vmax=3.)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=rvb)
    scalarMap._A=[]
    ax.view_init(azim=az,elev=el)
    ax.dist*=1.06
    s = ax.scatter3D(Delta.T[5],Delta.T[6],Delta.T[7],marker='o', lw = 0,s=15,color=scalarMap.to_rgba(Delta.T[8]),alpha=1.)
    s.set_edgecolors = s.set_facecolors = lambda *args:None
    ax.set_xlim(0,1100)
    ax.set_ylim(0,1000)
    ax.set_zlim(0,1000)
    ax.set_xlabel(r"$J_\lambda/{\rm\, kpc\,km\,s}^{-1}$")
    ax.set_ylabel(r"$J_\mu/{\rm\, kpc\,km\,s}^{-1}$")
    ax.set_zlabel(r"$J_\nu/{\rm\, kpc\,km\,s}^{-1}$")
    ax.xaxis._axinfo['ticklabel']['space_factor'] =1.1
    ax.xaxis._axinfo['label']['space_factor'] = 2.4
    ax.yaxis._axinfo['ticklabel']['space_factor'] =1.1
    ax.yaxis._axinfo['label']['space_factor'] = 2.4
    ax.zaxis._axinfo['ticklabel']['space_factor'] =1.
    ax.zaxis._axinfo['label']['space_factor'] = 2.4
    cbar_ax = f.add_axes([0.46 / 3., 0.87, 0.42 / 3., 0.03])
    ccc = f.colorbar(scalarMap,ticks=[0, 1, 2, 3],cax=cbar_ax,orientation='horizontal')
    cbar_ax.set_xticklabels(['Box','Short-axis\nloop', 'Inner \n long-axis loop', 'Outer long-axis\nloop'],fontsize=12)
    ccc.set_label(r'${\bf Orbit\,class}$',fontsize=12)
    ccc.ax.xaxis.set_ticks_position('top')
    ccc.ax.xaxis.set_label_position('top')

    rvb = make_colormap([c('black'), c(colA), 0.1, c(colA), c(colB), 0.5, c(colB),c(colC)])
    ax = f.add_subplot(1, 4, 2, projection='3d')
    color_max = np.max(Delta.T[2])
    color_min = np.min(Delta.T[2])
    # rvb = plt.get_cmap('brg')
    cNorm = colors.Normalize(vmin=color_min,vmax=color_max)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=rvb)
    scalarMap._A=[]
    ax.view_init(azim=az,elev=el)
    ax.dist*=1.06
    s = ax.scatter3D(Delta.T[5],Delta.T[6],Delta.T[7],marker='o', lw = 0,s=15,color=scalarMap.to_rgba(Delta.T[2]),alpha=1.)
    s.set_edgecolors = s.set_facecolors = lambda *args:None

    # grid_x, grid_y = np.mgrid[0:1100:100j, 0:1000:100j]
    # grid_z = griddata(Delta.T[5:7].T, Delta.T[7], (grid_x, grid_y), method='linear')
    # color_grid = griddata(Delta.T[5:7].T, Delta.T[3], (grid_x, grid_y), method='linear')
    # ax.plot_surface(
    #     grid_x, grid_y, grid_z, rstride=1, cstride=1,
    #     facecolors=scalarMap.to_rgba(color_grid),
    #     linewidth=0, antialiased=False, shade=False)

    ax.set_xlim(0,1100)
    ax.set_ylim(0,1000)
    ax.set_zlim(0,1000)
    ax.set_xlabel(r"$J_\lambda/{\rm\, kpc\,km\,s}^{-1}$")
    ax.set_ylabel(r"$J_\mu/{\rm\, kpc\,km\,s}^{-1}$")
    ax.set_zlabel(r"$J_\nu/{\rm\, kpc\,km\,s}^{-1}$")
    ax.xaxis._axinfo['ticklabel']['space_factor'] =1.1
    ax.xaxis._axinfo['label']['space_factor'] = 2.4
    ax.yaxis._axinfo['ticklabel']['space_factor'] =1.1
    ax.yaxis._axinfo['label']['space_factor'] = 2.4
    ax.zaxis._axinfo['ticklabel']['space_factor'] =1.
    ax.zaxis._axinfo['label']['space_factor'] = 2.4
    cbar_ax = f.add_axes([1.05 / 3., 0.87, 0.4 / 3., 0.03])
    ccc = f.colorbar(scalarMap,cax=cbar_ax,orientation='horizontal')
    ccc.set_label(r'$\Delta J_\lambda/{\rm\, kpc\,km\,s}^{-1}$')
    ccc.ax.xaxis.set_ticks_position('top')
    ccc.ax.xaxis.set_label_position('top')

    ax = f.add_subplot(1, 4, 3, projection='3d')
    color_max = np.max(Delta.T[3])
    color_min = np.min(Delta.T[3])
    cNorm = colors.Normalize(vmin=color_min,vmax=color_max)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=rvb)
    scalarMap._A=[]
    ax.view_init(azim=az,elev=el)
    ax.dist*=1.06
    s = ax.scatter3D(Delta.T[5],Delta.T[6],Delta.T[7],marker='o', lw = 0,s=15,color=scalarMap.to_rgba(Delta.T[3]),alpha=1.)
    s.set_edgecolors = s.set_facecolors = lambda *args:None
    ax.set_xlim(0,1100)
    ax.set_ylim(0,1000)
    ax.set_zlim(0,1000)
    ax.set_xlabel(r"$J_\lambda/{\rm\, kpc\,km\,s}^{-1}$")
    ax.set_ylabel(r"$J_\mu/{\rm\, kpc\,km\,s}^{-1}$")
    ax.set_zlabel(r"$J_\nu/{\rm\, kpc\,km\,s}^{-1}$")
    ax.xaxis._axinfo['ticklabel']['space_factor'] =1.1
    ax.xaxis._axinfo['label']['space_factor'] = 2.4
    ax.yaxis._axinfo['ticklabel']['space_factor'] =1.1
    ax.yaxis._axinfo['label']['space_factor'] = 2.4
    ax.zaxis._axinfo['ticklabel']['space_factor'] =1.
    ax.zaxis._axinfo['label']['space_factor'] = 2.4
    cbar_ax = f.add_axes([1-1.37/3., 0.87, 0.4 / 3., 0.03])
    ccc = f.colorbar(scalarMap,cax=cbar_ax,orientation='horizontal')
    ccc.set_label(r'$\Delta J_\mu/{\rm\, kpc\,km\,s}^{-1}$')
    ccc.ax.xaxis.set_ticks_position('top')
    ccc.ax.xaxis.set_label_position('top')

    ax = f.add_subplot(1, 4, 4, projection='3d')
    color_max = np.max(Delta.T[4])
    color_min = np.min(Delta.T[4])
    cNorm = colors.Normalize(vmin=color_min,vmax=color_max)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=rvb)
    scalarMap._A=[]
    ax.view_init(azim=az,elev=el)
    ax.dist*=1.06
    s = ax.scatter3D(Delta.T[5],Delta.T[6],Delta.T[7],marker='o', lw = 0,s=15,color=scalarMap.to_rgba(Delta.T[4]),alpha=1.)
    s.set_edgecolors = s.set_facecolors = lambda *args:None
    ax.set_xlim(0,1100)
    ax.set_ylim(0,1000)
    ax.set_zlim(0,1000)
    ax.set_xlabel(r"$J_\lambda/{\rm\, kpc\,km\,s}^{-1}$")
    ax.set_ylabel(r"$J_\mu/{\rm\, kpc\,km\,s}^{-1}$")
    ax.set_zlabel(r"$J_\nu/{\rm\, kpc\,km\,s}^{-1}$")
    ax.xaxis._axinfo['ticklabel']['space_factor'] =1.1
    ax.xaxis._axinfo['label']['space_factor'] = 2.4
    ax.yaxis._axinfo['ticklabel']['space_factor'] =1.1
    ax.yaxis._axinfo['label']['space_factor'] = 2.4
    ax.zaxis._axinfo['ticklabel']['space_factor'] =1.
    ax.zaxis._axinfo['label']['space_factor'] = 2.4
    cbar_ax = f.add_axes([1 - 0.78 / 3., 0.87, 0.4 / 3., 0.03])
    ccc = f.colorbar(scalarMap,cax=cbar_ax,orientation='horizontal')
    ccc.set_label(r'$\Delta J_\nu/{\rm\, kpc\,km\,s}^{-1}$')
    ccc.ax.xaxis.set_ticks_position('top')
    ccc.ax.xaxis.set_label_position('top')

    plt.savefig('3d_error_diagram.pdf',bbox_inches='tight')
    plt.clf()

# plot_closed_short()
# plot_closed_long()
# Delta_withEnergy()
DeltaCheck()
orbit_accuracy_log10fractional()
orbit_accuracy_sum()
orbit_accuracy_absolute()
# testLoop("testloop_v1.dat",np.array([54.704069797,752.36399863,77.8328506573]),'testLoop.png',2,1,'Short-axis loop')
# testLoop("testloop_v2.dat",np.array([704.121354985,176.648233384,139.311624494]),'testLoop_v2.png',1,0,'Box')
# testLoop("testloop_v3.dat",np.array([49.7338329784,103.389963613,678.943849563]),'testLoop_v3.png',2,0,'Long-axis loop')

# plot_SoS("testloop_v1_sos.dat","testloop_v1_sos2.dat",'sos_v1.pdf',r'$x$',r'$_x$','Short-axis loop, ' + r'$y=0,\,z=0$')
# # plot_SoS("testloop_v1_sos_y.dat","testloop_v1_sos2_y.dat",'sos_v1_y.pdf',r'$y$',r'$_y$','Short-axis loop, ' + r'$x=0,\,z=0$')
# plot_SoS("testloop_v2_sos.dat","testloop_v2_sos2.dat",'sos_v2.pdf',r'$x$',r'$_x$','Box, ' + r'$y=0,\,z=0$')
# plot_SoS("testloop_v3_sos.dat","testloop_v3_sos2.dat",'sos_v3.pdf',r'$y$',r'$_y$','Long-axis loop, ' + r'$x=0,\,z=0$')

# plot_jeans_solution('jeans_a3.28_0.dat','jeans_a328_0.pdf',[8],[0],r'$x=x_i,y=1{\,\rm kpc},z=1{\,\rm kpc}$', r'Radial bias, $\zeta=3.28$')
# plot_jeans_solution('jeans_a3.28_1.dat','jeans_a328_1.pdf',[9],[1],r'$x=1{\,\rm kpc},y=x_i,z=1{\,\rm kpc}$', '')
# plot_jeans_solution('jeans_a3.28_2.dat','jeans_a328_2.pdf',[10],[2],r'$x=1{\,\rm kpc},y=1{\,\rm kpc},z=x_i$', '')
# plot_jeans_solution('jeans_a3.28_3.dat','jeans_a328_3.pdf',[8],[0],r'$x=x_i,y=x/2,z=x/3$', r'Radial bias, $\zeta=3.28$')
# plot_jeans_solution('jeans_a3.28_3.dat','jeans_a328_4.pdf',[9],[1],r'$x=x_i,y=x/2,z=x/3$', '')
# plot_jeans_solution('jeans_a3.28_3.dat','jeans_a328_5.pdf',[10],[2],r'$x=x_i,y=x/2,z=x/3$', '')

# plot_jeans_solution('jeans_a0.7_0.dat','jeans_a07_0.pdf',[8],[0],r'$x=x_i,y=1{\,\rm kpc},z=1{\,\rm kpc}$', r'Tangential bias, $\zeta=0.7$')
# plot_jeans_solution('jeans_a0.7_1.dat','jeans_a07_1.pdf',[9],[1],r'$x=1{\,\rm kpc},y=x_i,z=1{\,\rm kpc}$', '')
# plot_jeans_solution('jeans_a0.7_2.dat','jeans_a07_2.pdf',[10],[2],r'$x=1{\,\rm kpc},y=1{\,\rm kpc},z=x_i$', '')
# plot_jeans_solution('jeans_a0.7_3.dat','jeans_a07_3.pdf',[8],[0],r'$x=x_i,y=x/2,z=x/3$', r'Tangential bias, $\zeta=0.7$')
# plot_jeans_solution('jeans_a0.7_3.dat','jeans_a07_4.pdf',[9],[1],r'$x=x_i,y=x/2,z=x/3$', '')
# plot_jeans_solution('jeans_a0.7_3.dat','jeans_a07_5.pdf',[10],[2],r'$x=x_i,y=x/2,z=x/3$', '')

# twoDdensity('tangential_bias_density.dat','a0.700000_density_split.dat',r'$x=x_i,y=1{\,\rm kpc},z=1{\,\rm kpc}$',r'Tangential bias, $\zeta=0.7$','density_tangential.pdf')
# twoDdensity('radial_bias_density.dat','a3.280000_density_split.dat',r'$x=x_i,y=1{\,\rm kpc},z=1{\,\rm kpc}$',r'Radial bias, $\zeta=3.28$','density_radial.pdf')

# error_planes()

# constituent_density('density_1.dat','density_1.eps',r'$x=1{\,\rm kpc},y=x_i,z=1{\,\rm kpc}$')
# constituent_density('density_2.dat','density_2.eps',r'$x=1{\,\rm kpc},y=1{\,\rm kpc},z=x_i$')
# constituent_density('density_3.dat','density_3.eps',r'$x=x_i,y=x/2,z=x/3$')
