#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
from scipy import sparse
import sys
import time
from scipy.stats.mstats import mquantiles
from scipy.interpolate import interp1d,interp2d,RectBivariateSpline
import os,sys
#sys.path.append('/afs/ipp/home/t/todstrci/TRANSP/')

#import equ_map 

from scipy import constants
#sys.path.append('/afs/ipp/home/g/git/python/repository/')
#sys.path.append('/afs/ipp/u/transp/tr_client/SUN/test/')

#import ufiles as uf
#from plot_profile import PlotProfileAnim,PlotProfile
#kk = equ_map.equ_map()



def LoadAtomData( file):
    f = open(file)
    header = f.readline()
    n_ions, n_ne, n_te = header.split()[:3]
    details = ' '.join(header.split()[3:])

    f.readline() #remove splitter

    n_ions,n_ne,n_te = int(n_ions),int(n_ne),int(n_te)
    Te = []
    ne = []
    while len(ne)< n_ne:
        line = f.readline()
        ne = ne+[float(n) for n in line.split()]
    while len(Te)< n_te:
        line = f.readline()
        Te = Te+[float(t) for t in line.split()]  
    Te = array(Te)
    ne = array(ne)
    ions_data = []
    for i_ion in range(n_ions):
        f.readline()#remove splitter
        plsx = []
        while len(plsx)< n_ne*n_te:
            line = f.readline()
            plsx = plsx+[float(L) for L in line.split()] 
        plsx = array(plsx).reshape(n_te, n_ne)
        ions_data.append(array(plsx))
    
    return ne,Te,ions_data


def LoadFile( file):
    f = open(file)
    header = f.readline()
    n_te = header.split()[0]

    Te_label = f.readline().split()[1] #remove splitter
    fact = 1e3 if Te_label=='[keV]' else 1 
    n_te = int(n_te)
    Te = []
    while len(Te)< n_te:
        line = f.readline()
        Te = Te+[float(t) for t in line.split()]  
    Te = array(Te)
    ions_data = {}
    head = f.readline()

    while head != '':
        plsx = []
        while len(plsx)< n_te:
            line = f.readline()
            plsx.extend([float(L) for L in line.split()])
        element = head[:2].strip(' _')
        ions_data[element] = array(plsx) *fact**2/1e6
        head = f.readline()
        
    
    return Te*fact,ions_data




def CalcBackgroundRadiation(ne_data,te_data,Zeff,filt=5,main_ion='D',version='puetti_'):
    #assuming only D and B
    #lz_sxr_5  75um Be filter
    #lz_sxr_6  75um Be filter old data!!
    #lz_sxr_7  20um Be filter
    #lz_sxr_8  0um Be filter
    #lz_sxr_10  00um AXUV

    #lz  no  Be filter

    path = os.path.dirname(os.path.realpath(__file__))
    
    file = 'lz_%sbolo.dat'%version if filt is None else 'lz_%ssxr_%d.dat'%(version,filt)   

    te,coeff = LoadFile(os.path.join(path,file))
    impurity = 'B'
    Z_i = CalcZ(te_data,element='B')

    if not impurity in coeff:
        Z_i = 6
        impurity = 'C'
    if main_ion in ['D','H']:    
        Z_m = 1.
        main_ion =  'H'
    if main_ion == 'He':
        Z_m = 2.

    
    c_m = abs((Z_i-Zeff)/(Z_m**2-Z_m*Z_i))
    c_i = abs((Zeff-Z_m)/(Z_i**2-Z_i*Z_m))
    #print 'filt', type(filt)
    

    te_data = maximum(te_data,50)  #troubles with values close to zero
    ne_data = maximum(ne_data,50)  #troubles with values close to zero

    lte = log10(te_data)
    lne = log10(ne_data)


    #te[te]
    LzM = interp(lte, log10(te),log10(coeff[main_ion]))

    emiss =  10**(LzM+2*lne)*c_m  #W/m^3



    Lzi = interp(lte, log10(te),log10(coeff[impurity]))

    emiss +=  10**(Lzi+2*lne)*c_i  #W/m^3
    
    #if filt == 10: emiss/= 0.27#A/W
    
    #te_corr, corr = loadtxt(os.path.join(path,'Bremsstrahlung_Lz_correction.txt')).T
    #print 'bremsstralung correction was applied!'
    
    #emiss *= interp(te_data, te_corr,corr)
    if filt in [5,6,7]:
        emiss *= 2


    return emiss





    #ne_data=ne_data/1e6  # m^-3 => cm^-3

    
    #file = 'prsx5_h.dat'  #hydrogen data, applied Be filter for SXR, recombination bremstrahlung
    #ne,te,coeff = LoadAtomData( file)  #cm^-3,eV, ? 
    #int_fun = RectBivariateSpline(te,ne,coeff[0],kx=3, ky=3, s=0)
    
    #prsx5 = int_fun.ev( log10(te_data).ravel(),log10(ne_data).ravel())
    #nr,nt = te_data.shape
    #prsx5 = prsx5.reshape(nr,nt)  #cooling coefficient
    
    #emiss =  10**prsx5*ne_data**2*c_D*1e6  #W/m^3
    ##plot(prsx5.T,'b')

    ##plot(te,coeff[0],'k',label='h r')


    #file = 'prsx5_b.dat'  #boron data, applied Be filter for SXR, recombination bremstrahlung

    #ne,te,coeff = LoadAtomData( file)
    
    #int_fun = RectBivariateSpline(te,ne,coeff[-1],kx=3, ky=3, s=0)  #BUG!!! je tu jen poslední iont! 
    
    #prsx5 = int_fun.ev( log10(te_data).ravel(),log10(ne_data).ravel())
    
    #nr,nt = te_data.shape
    #prsx5 = prsx5.reshape(nr,nt)  #cooling coefficient
   

    #emiss +=  10**prsx5*ne_data**2*c_B*1e6  #W/m^3
    
    #emiss2 = 10**prsx5*ne_data**2*c_B*1e6
    #plot(prsx5.T,'r')

    
    #plot((10**(prsx5)*ne_data**2*c_B*1e6).T ,'r')
    #show()
    #exit()
    #plot(te,coeff[-4],'g',label='b r')

    #plot(te,coeff[-3],'g',label='b r')

    #plot(te,coeff[-2],'g',label='b r')

    #plot(te,coeff[-1],'g',label='b r')

 
    return emiss
    
    
    
def CalcZ(te_data,element='W'):
    #assuming only W
    
    
    #import IPython
    #IPython.embed()
    path = os.path.dirname(os.path.realpath(__file__))

    file = 'lz_puetti_meanz.dat'   

    te,coeff = LoadFile(os.path.join(path,file))
    te_data = maximum(te_data,1)


    Z = coeff[element] *1e6  #1e6 is a  "bug" in the  LoadFile
     
    Z=interp(log10(te_data), log10(te),Z)  
    
    
    return Z
    
def CalcRadiation(ne_data,te_data,c_imp,filt=5,element='W',version=None):
    #assuming only W
    #lz_sxr_5  75um Be filter
    #lz_sxr_7  20um Be filter
    #lz_sxr_6  75um Be filter old data!!
    #lz_sxr_8  75um Be origin L_W 
    #version '' - old
    #version 'puetti_' - new by Thomas
    #version None - by Thomas+improved SXR Lz
    #filt 10 - axuv 
    path = os.path.dirname(os.path.realpath(__file__))

    #file = 'lz.dat' if filt is None else 'lz_sxr_%d.dat'%filt  
    load_new=False
    if version is None:
        load_new = True
        version = 'puetti_'
    file = 'lz_%sbolo.dat'%version if filt is None else 'lz_%ssxr_%d.dat'%(version,filt)   


    te,coeff = LoadFile(os.path.join(path,file))

    #experimental values of the cooling coefficient
    if element == 'W' and filt == 5 and load_new:
        #print 'new Lz coeff!!'
        te, coeff = loadtxt(os.path.join(path,'new_W_lz2.txt'),unpack=True)
        coeff = {'W':coeff }

    te_data = maximum(te_data,te.min())
    Lz=interp(log10(te_data), log10(te),log10(coeff[element]))
    emiss = 10**(Lz+2*log10(ne_data))*c_imp #W/m^3

    if filt == 10: emiss/= 0.27#A/W

    return emiss
    
    
    
    
    
 
    
def plot_cooling(elements,file_name):
    
    
    import fconf
    fconf.useslide()
    import matplotlib.pylab as plt
    file = 'lz_sxr_5.dat' 
    #elements = ['W','Ar','Ni','Si']
    color=['k-','r--','b-.','g:','y','m']

    te,coeff = LoadFile(file)
    #print coeff.keys()
    plt.figure(figsize=(5,5))
    plt.clf()
    
    for c,name in zip(color,elements):
        plt.loglog(te*1e3,coeff[name]/amax(coeff[name])*1e-30,c,label=name )
        
    #file = 'lz_sxr_7.dat' 
    #te,coeff = LoadFile(file)
    #for c,name in  zip(color,elements):
        #plt.loglog(te,coeff[name],c+'--')
    
    plt.xlim(100,1e4)
    plt.ylim(1e-35,None)
    plt.ylabel('$L_z $ [W m$^3$]')
    plt.xlabel('$T_e$ [eV]')

    plt.legend(loc='lower right')

    plt.savefig(file_name+'.png')
    plt.savefig(file_name+'.pdf',bbox_inches='tight' )
    plt.savefig(file_name+'.svg',bbox_inches='tight' )

    plt.show()
    
    
    
