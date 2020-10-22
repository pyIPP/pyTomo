import numpy as np
import MDSplus as mds
from scipy.interpolate import interp1d
from scipy.io import loadmat
from matplotlib.pylab import *


def calc_zeff(shot,Z=[6,8],pZ=[.8,.2]):
    # Estimate the Zeff profile from Thomson and XTOMO data.
    # XTOMO data are inverted with GTI package on a Flux1D grid.
    # INPUT:    - shot: shot number. From shot, the operating gas is determined
    #           and therefore the charge of the background ions b.
    #           - Z: the charge of the impurity ions. default: [6 8] (carbon
    #           and oxygen)
    #           - pZ: the concentration of impurities. default [0.8 0.2].
    #
    # OUTPUT: out is a structure with fields   
    #           - Zx: effective charge profile [time x rho].
    #           - nimp: overall impurity density profile [time x rho]
    #           - cimp: overall impurity concentration profile [time x rho]
    #           - t: time (Thomson times actually)
    #           - rho: normalized poloidal flux
    #           - ne: density profile from Thomson
    #           - Te: temperature profile from Thomson
    #
    # WARNING:  - Impurity injection is not implemented in this version.
    #           - The method is not reliable for Te<200 eV.
    #           - Only 47 mu Be filter taken into account in this version.
    #
    # CRPP/EPFL  April 2012   Benoit Labit



    if np.size(pZ)!=np.size(Z):
        raise Exception('Z and pZ must have the same size')

    # Operating gas
    if shot>=42146 and shot<=42227:
        gas_type = 'H2'; # gas type
        b = 1;
    elif shot>=42495 and shot<=42524:
        gas_type = 'He';
        b=2;
    else:
        gas_type='D2'; # Default gas type
        b=1;
    
    Z = asarray(Z)
    
    #BUG!!!
    mds_server = 'localhost:8001'
    c= mds.Connection(mds_server )
    c.openTree('tcv_shot',shot)


    c.get('$shot');
    # Download density and temperature profiles from Thomson
    print(' download Thomson profiles')
    ne0 = np.asarray(c.get(r'\results::thomson.profiles.auto:ne'));
    Te0 = np.asarray(c.get(r'\results::thomson.profiles.auto:te'));
    rhoT = np.asarray(c.get(r'\results::thomson:profiles:auto:rho')); #			-	TCV standard poloidal rho=sqrt(psi_bar)= linspace(0,1,41)
    tT = np.asarray(c.get(r'\results::thomson.profiles.auto:time'));# s	time vector (usually shorter than time vector of raw data) ne.dim{1};
    #print ne0
    ind = ~np.all(np.isnan(ne0),1)
    ne0 = ne0[ind]
    Te0 = Te0[ind]
    rhoT = rhoT[ind]
    

    c.closeTree('tcv_shot',shot)

   
    # Tomography inversion 
    #PWD = pwd;
    print('.  tomography for XTOMO')
    #gti = xtomo_gti(shot,tT,[],3);
    #cd(PWD);
    #BUG

    #data = np.loadtxt('/mnt/hardisk/tomo5/tmp/emiss_profile_%d.txt'%shot);
    data = np.loadtxt('/mnt/harddisk/tomo5/tmp/emiss_profile_%d.txt'%shot)
    Emix0 = data[:,1:]
    tvec = data[:,0]
    
    ind = ~np.isnan(tT)&(tT>=tvec[0])&(tT<=tvec[-1])
    ne0 = ne0[:,ind]
    Te0 = Te0[:,ind]
    tT = tT[ind]
    
    
    rho = np.linspace(0,1,Emix0.shape[1])
    Emix0 = interp1d(tvec,Emix0,axis=0)(tT )
    Emix0 = interp1d(rho,Emix0,axis=1)(rhoT ).T



    # Compute normalized emissivity (based on IONEQ)

    thickness = 47; #Be filt thickness
    cimp=[]; Zx=[];

    emix = zxpro_emix(b,Te0,thickness);
    deadp = loadmat('zxpro_deadcorr.mat')['deadp'].flatten()
    dead = np.polyval(deadp,np.log(Te0));
    emix *= dead;
    D=zeros_like(emix)
    for z,pz in zip(Z,pZ):
        E = dead*zxpro_emix(z,Te0,thickness)
        D+=pz*(E-z*emix)
    
    import IPython
    IPython.embed()
    
    plot( Emix0.mean(1))
    plot(np.nanmean( ne0**2*emix,1))
    show()
    
    
    
    cimp = (b*Emix0/(ne0**2)-emix)/D; # Eq. 4.27 Ivo's thesis
    Zx = b + cimp*sum(Z*(Z-b)*pZ); # Eq. 4.28 Ivo's thesis

    # The method is not reliable for low temperature
    TeMin = 200;
    Zx[Te0<TeMin]=nan;
    cimp[Te0<TeMin]=nan;
    cimp[Zx>inner(Z,pZ)]=nan;
    Zx[Zx>inner(Z,pZ)]=nan;
    # Outputs
    class out_:pass
    out = out_()
    out.Zx=Zx; out.rho = rhoT; out.nimp = cimp*ne0;
    out.t=tT; out.ne = ne0; out.Te = Te0; out.cimp = cimp;

    return out


def zxpro_emix(Z,te,thick=47):

    ####################################################################
    # function emix=emix(Z,te,Bethick)
    #   supplies normalized emissivity for z and given te vector
    #   preferred use with global definition in main program to avoid 
    #   frequent disk access:
    #         global ION_Z POLYORDER_0 POLYORDER_47 IONPOLY_0 IONPOLY_47
    #        load ionmatrix
    #####################################################################

    scale=4e-32; # scaling from IONEQ defaults
    #global ION_Z POLYORDER_0 POLYORDER_47 IONPOLY_0 IONPOLY_47 

    #if len(ION_Z)==0:
    data = loadmat('zxpro_ionmatrix') # casual use, normally not necessary
    #print data.keys()
    index= np.where(data['ION_Z']==Z)[0][0];
    if thick!=47  and thick!=0:
        print(('EMIX error: requested filter thickness unavailable'+str(thick)+' available:0, 47')) 
    
    ok = (te>=200) & (te<=20000) & np.isfinite(te)
    emix=np.nan*np.ones_like(te);
    #import IPython
    #IPython.embed()
    if thick==47:
        p=data['IONPOLY_47'][index, -int(data['POLYORDER_47'][0,index])-1:]
    else:
        p=data['IONPOLY_0'][index, -int(data['POLYORDER_0'][0,index])-1:]
    

    emix[ok]=scale*np.exp(np.polyval(p,np.log(te[ok])))
    
    return emix





if __name__ == '__main__':
    
    
    out = calc_zeff(40412)
    import matplotlib.pylab as plt
    plt.plot(out.t, out.Zx)
    plt.show()
    

