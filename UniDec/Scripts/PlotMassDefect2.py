__author__ = 'Michael.Marty'
__author__ = 'Michael.Marty'

import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import matplotlib as mpl
from matplotlib import rcParams

from unidec_modules.isolated_packages import option_d

rcParams['ps.useafm']=True
rcParams['ps.fonttype']=42
rcParams['pdf.fonttype']=42
rcParams['lines.linewidth']=0.5
rcParams['axes.linewidth']=0.5
#rcParams['font.size']=18

def mergedata(data1,data2):
    f = interp1d(data2[:, 0], data2[:, 1],bounds_error=False,fill_value=0)
    intx=data1[:,0]
    inty=f(intx)
    newdat = np.column_stack((intx, inty))
    return newdat

def mergedata2d(X1,Y1,X2,Y2,Z2):
    #f=interp2d(np.unique(np.ravel(X2)),np.unique(np.ravel(Y2)),Z2,bounds_error=False,fill_value=0)
    f=interp2d(np.ravel(X2),np.ravel(Y2),np.ravel(Z2),bounds_error=False,fill_value=0)
    newx=np.unique(np.ravel(X1))
    newy=np.unique(np.ravel(Y1))
    Zout=f(newx,newy)
    return Zout

topdir="C:\\Data\\"
proteins=["AmtB","AqpZ"]
pmass=[126783,98848]
lmass=[760.076,760.076]

datafiles2=[["AmtB_POPC\\AmtB_04","AmtB_POPC\\AmtB_10_E3","AmtB_DMPC\\AmtB_07"],["AqpZ_POPC\\Aqpz_07","AqpZ_POPC\\Aqpz_11_E3","AqpZ_DMPC\\Aqpz_09_2"]]
datafiles=[[["AmtB_POPC\\AmtB_02","AmtB_POPC\\AmtB_04","AmtB_POPC\\AmtB_05"]
               ,["AmtB_POPC\\AmtB_10_E3","AmtB_POPC\\AmtB_12_E3","AmtB_POPC\\AmtB_13_E3"]
               ,["AmtB_DMPC\\AmtB_07","AmtB_DMPC\\AmtB_08_2","AmtB_DMPC\\AmtB_09_2"]
            ],
           [["AqpZ_POPC\\Aqpz_05_2","AqpZ_POPC\\Aqpz_06_2","AqpZ_POPC\\Aqpz_07"]
               ,["AqpZ_POPC\\Aqpz_11_E3","AqpZ_POPC\\Aqpz_12_E3","AqpZ_POPC\\Aqpz_13_E3"]
               ,["AqpZ_DMPC\\Aqpz_08","AqpZ_DMPC\\Aqpz_09_2","AqpZ_DMPC\\Aqpz_10"]
            ]]
ncols=len(datafiles2[0])+1
ratios=np.ones(ncols)
ratios[-1]=0.05

gs=gridspec.GridSpec(2,ncols,width_ratios=ratios)
gs.update(wspace=0.1,hspace=0.1)

fig=plt.figure(figsize=(6.3,4))

#magicnumber=[[12.882,20.045,72.040,54,10.32],[7.498,16.385,78.802,44,5.58]]
#magicnumber=[[10.32,20.045,72.040,54],[5.58,16.385,78.802,44]]
magicnumber=[[7.498,78.002],[12.882,72.040]]

cmap= option_d.viridis
#cmap="jet"

size=8
for i,protein in enumerate(proteins):

    maximum=220
    minorticks=np.arange(0,maximum,25)
    majorticks=np.arange(0,1,0.2)
    #MS PLOTS

    group=datafiles[i]
    for j,group2 in enumerate(group):

        #plt.subplot(2,2,i+1)
        plt.subplot(gs[i*ncols+j])


        ax=plt.gca()
        ax.set_xticks(minorticks,minor=True)
        #ax.set_yticks(majorticks)
        ax.set_xlim(20,maximum)
        ax.set_ylim(0,1)

        if i==0:
            ax.set_xticklabels([])
        if j>0:
            ax.set_yticklabels([])
        else:
            #ax.set_yticklabels([0,"",1],fontsize=size)
            if i==0:
                plt.ylabel("AmtB\nNormalized Mass Defect",fontsize=size)
            else:
                plt.ylabel("AqpZ\nNormalized Mass Defect",fontsize=size)
        plt.setp(ax.get_yticklabels(),fontsize=size)
        plt.setp(ax.get_xticklabels(),fontsize=size,color="k")
        labels=["POPC","POPC MSP1E3D1(-)","DMPC","DMPC E3"]
        styles=["-",":","--","-."]
        if i==0:
            plt.title(labels[j],fontsize=size)
        datall=[]
        switch=0
        for folder in group2:
            path=os.path.join(topdir,os.path.join(folder,"Total_2D_Mass_Defects.txt"))
            print path
            data=np.loadtxt(path)
            data[:,2]=data[:,2]/np.amax(data[:,2])

            X = data[:, 0]/1000.
            Y = data[:, 1]
            C = data[:, 2]

            #bool1=X>100
            #C=C/np.amax(C[bool1])
            #C=np.clip(C,0,1)


            xvals = np.unique(X)
            yvals = np.unique(Y)
            xlen = len(xvals)
            ylen = len(yvals)
            newgrid = np.reshape(C, (xlen, ylen))
            print newgrid.shape
            newx = np.reshape(X, (xlen, ylen))
            newy = np.reshape(Y, (xlen, ylen))
            datall.append(newgrid)
            '''
            if switch==0:
                datall.append(newgrid)
                switch=1
                topx=newx
                topy=newy
            else:
                datall.append(mergedata2d(topx,topy,newx,newy,newgrid))
            '''

        newgrid=np.sum(datall,axis=0)
        print newgrid.shape
        cax=plt.contourf(newx, newy, newgrid, 100, cmap=cmap)
        ticcol="w"
        ax=plt.gca()
        for line in ax.xaxis.get_ticklines():
            line.set_color(ticcol)
        for line in ax.yaxis.get_ticklines():
            line.set_color(ticcol)
        #ax.tick_params(axis="x",which="both",colors=ticcol)
        if i>0:
            plt.xlabel("Mass (kDa)",fontsize=size)
        for v in magicnumber[i]:
            pass
            #plt.vlines(v,0,1)

    plt.subplot(gs[i*ncols+ncols-1])
    ax=plt.gca()
    cbar=mpl.colorbar.ColorbarBase(ax,cmap=cmap,orientation="vertical",ticks=[0,0.5,1])
    #cbar=plt.colorbar(use_gridspec=False,ticks=[0,np.amax(newgrid)/2,np.amax(newgrid)])
    cbar.ax.get_yaxis().set_tick_params(direction='out')
    cbar.ax.set_yticklabels(["0","%","100"],fontsize=size)


plt.savefig("C:\\Users\\michael.marty\\Dropbox\\Research\\MSNanodisc\\Figures\\Figure_All_Kendrick.png")
plt.savefig("C:\\Users\\michael.marty\\Dropbox\\Research\\MSNanodisc\\Figures\\Figure_All_Kendrick.pdf")
plt.show()

