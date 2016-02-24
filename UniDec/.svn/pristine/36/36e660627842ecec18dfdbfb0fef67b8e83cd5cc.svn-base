__author__ = 'Michael.Marty'


__author__ = 'Michael.Marty'

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
from matplotlib import rcParams


rcParams['ps.useafm']=True
rcParams['ps.fonttype']=42
rcParams['pdf.fonttype']=42
rcParams['lines.linewidth']=0.5
rcParams['axes.linewidth']=0.5
rcParams['font.size']=8

def mergedata(data1,data2):
    f = interp1d(data2[:, 0], data2[:, 1],bounds_error=False,fill_value=0)
    intx=data1[:,0]
    inty=f(intx)
    newdat = np.column_stack((intx, inty))
    return newdat

topdir="C:\\Data\\"
proteins=["AmtB","AqpZ"]

datafiles=[[["AmtB_POPC\\AmtB_02","AmtB_POPC\\AmtB_04","AmtB_POPC\\AmtB_05"],["AmtB_POPC\\AmtB_10_E3","AmtB_POPC\\AmtB_12_E3","AmtB_POPC\\AmtB_13_E3"],["AmtB_DMPC\\AmtB_07","AmtB_DMPC\\AmtB_08_2","AmtB_DMPC\\AmtB_09_2"]],
           [["AqpZ_POPC\\Aqpz_05_2","AqpZ_POPC\\Aqpz_06_2","AqpZ_POPC\\Aqpz_07"],["AqpZ_POPC\\Aqpz_11_E3","AqpZ_POPC\\Aqpz_12_E3","AqpZ_POPC\\Aqpz_13_E3"],["AqpZ_DMPC\\Aqpz_08","AqpZ_DMPC\\Aqpz_09_2","AqpZ_DMPC\\Aqpz_10"]]]
    #,["AmtB_DMPC\\AmtB_11_E3"]]]


gs=gridspec.GridSpec(2,2)
gs.update(wspace=0.1,hspace=0,right=0.85)

fig=plt.figure(figsize=(4.13,3))

#magicnumber=[[12.882,20.045,72.040,54,10.32],[7.498,16.385,78.802,44,5.58]]
#magicnumber=[[10.32,20.045,72.040,54],[5.58,16.385,78.802,44]]
magicnumber=[[7.498,78.002],[12.882,72.040]]

colors=["silver","gray","k"]
size=8
dattot=[]
for i,protein in enumerate(proteins):

    maximum=120
    #plt.subplot(2,2,i+3)
    plt.subplot(gs[i+2])
    ax=plt.gca()
    ax.set_xlim(0,maximum)
    minorticks=np.arange(0,maximum,5)
    majorticks=np.arange(0,1.2,0.2)
    ax.set_xticks(np.arange(0,120,20))
    ax.set_xticks(minorticks,minor=True)
    ax.set_yticks(majorticks)
    plt.xlabel("# of Lipids",fontsize=size)
    plt.setp(ax.get_xticklabels(),fontsize=size)
    if i==1:
        ax.set_yticklabels([])
    else:
        ax.set_yticklabels(["",0.2,0.4,0.6,0.8,1],fontsize=size)
        plt.ylabel("MD Probability",fontsize=size)


    dir="C:\\VMD\\"+protein+"\\AnalysisAlignedLong"
    os.chdir(dir)
    histall=[]
    for c in xrange(0,6):
        outfile=protein+"_class_hist_"+str(c)+".txt"
        hist=np.loadtxt(outfile)
        hist[:,1]=hist[:,1]/np.amax(hist[:,1])
        hist[:,1]=hist[:,1]-0
        histall.append(hist)
        if c==3:
            plt.fill_between(hist[:,0],hist[:,1]*0,hist[:,1],color=colors[0],label="Ionic")
        if c==4:
            pass
            #plt.plot(hist[:,0],hist[:,1],color="k",label="Ionic+\nH-Bond")
        if c==5:
            plt.fill_between(hist[:,0],hist[:,1]*0,hist[:,1],color=colors[2],label="Annular\nBelt")

    dir="C:\\VMD\\"+protein+"\\AnalysisAligned"
    os.chdir(dir)
    histall=[]
    for c in xrange(0,6):
        outfile=protein+"_class_hist_head_"+str(c)+".txt"
        hist=np.loadtxt(outfile)
        hist[:,1]=hist[:,1]/np.amax(hist[:,1])
        hist[:,1]=hist[:,1]-0
        histall.append(hist)
        if c==3:
            pass
            #plt.fill_between(hist[:,0],hist[:,1]*0,hist[:,1],color="gray",label="Ionic")
        if c==4:
            pass
            #plt.plot(hist[:,0],hist[:,1],color="k",label="Ionic+\nH-Bond")
        if c==5:
            plt.fill_between(hist[:,0],hist[:,1]*0,hist[:,1],color=colors[1],label="Head\nGroups")

    if i==1:
        plt.plot([], [], color=colors[2], linewidth=3,label="Annular\nBelt")
        plt.plot([], [], color=colors[1], linewidth=3,label="Head\nGroups")
        plt.plot([], [], color=colors[0], linewidth=3,label="Ionic\nContacts")


        leg=plt.legend(bbox_to_anchor=(1.30,0.85),loc=1,borderaxespad=0,fontsize=size)
        leg.get_frame().set_linewidth(0.5)
    else:
        pass
        #plt.text(1,0.9,"MD",fontsize=16)

    #MS PLOTS

    #plt.subplot(2,2,i+1)
    plt.subplot(gs[i])
    plt.title(protein,fontsize=size)
    ax=plt.gca()
    ax.set_xticks(minorticks,minor=True)
    ax.set_yticks(majorticks)
    ax.set_xlim(0,maximum)
    ax.set_xticklabels([])
    if i==1:
        ax.set_yticklabels([])
    else:
        ax.set_yticklabels(["",0.2,0.4,0.6,0.8,1],fontsize=size)
        plt.ylabel("MS Intensity",fontsize=size)

    labels=["POPC","POPC E3","DMPC","DMPC E3"]
    styles=["-",":","-",":"]
    c2=["b","b","r","r"]
    for k,dtype in enumerate(datafiles[i]):
        datall=[]
        maxlen=0
        for j,folder in enumerate(dtype):
            path=os.path.join(topdir,os.path.join(folder,"Extract_total_2D_Extract.txt"))
            print path
            data=np.loadtxt(path)
            data=data[data[:,0]<maximum]
            #data[:,1]=data[:,1]-np.amin(data[:,1])
            data[:,1]=data[:,1]/np.amax(data[:,1])
            datall.append(data)
            #plt.plot(data[:,0],data[:,1],color="b",linestyle=styles[j],label=str(j))
            print data.shape,folder
            if len(data)>maxlen:
                datafirst=data
                maxlen=len(data)
        print maxlen
        datall2=[]
        for data in datall:
            datall2.append(mergedata(datafirst,data))
        datall=np.array(datall2)
        datavg=np.average(datall[:,:,1],axis=0)
        datavg=datavg/np.amax(datavg)
        dattot.append(datavg)
        plt.plot(datafirst[:,0],datavg,color=c2[k],label=labels[k],linestyle=styles[k])



    for v in magicnumber[i]:
        pass
        #plt.vlines(v,0,1)

    if i==1:
        pass
        leg=plt.legend(bbox_to_anchor=(1.30,0.9),loc=1,borderaxespad=0,fontsize=size)#,title="Replication")
        leg.get_frame().set_linewidth(0.5)
    else:
        pass
        #plt.text(1,0.9,"MS",fontsize=16)
np.savetxt("C:\\Data\\lipidextracts.txt",dattot)
#plt.savefig("C:\\Users\\michael.marty\\Dropbox\\Research\\MSNanodisc\\Figures\\Figure_S10.png")
#plt.savefig("C:\\Users\\michael.marty\\Dropbox\\Research\\MSNanodisc\\Figures\\Figure_S10.pdf")

plt.show()