



import numpy as np
import MDSplus as mds
from scipy.interpolate import interp1d
from scipy.io import loadmat



shot =   42314
mds_server = 'localhost:8001'
c= mds.Connection(mds_server )
c.openTree('tcv_shot',shot)
print(c.get('$shot'));

time = c.get(r'\results::psi')

# Download density and temperature profiles from Thomson
print('.  download Thomson profiles')
Te2 =  c.get(r'\results::te_thomson')
Te = c.get(r'\results::thomson:te');



tvec = np.asarray(c.get(r'\results::thomson:times'));
Te = np.asarray(c.get(r'\results::thomson:te'));
Te_err = np.asarray(c.get(r'\results::thomson:te:error_bar'));
ne = np.asarray(c.get(r'\results::thomson:ne'));
ne_err = np.asarray(c.get(r'\results::thomson:ne:error_bar'));

psi = np.asarray(c.get(r'\results::thomson:psiscatvol'));
psimax = np.asarray(c.get(r'\results::thomson:psi_max'));
psi = np.sqrt(1-psi/psimax)

ne0 = np.asarray(c.get(r'\results::thomson.profiles.auto:ne'));
ne0_err = np.asarray(c.get(r'\results::thomson.profiles.auto:ne:error_bar'));
Te0 = np.asarray(c.get(r'\results::thomson.profiles.auto:te'));
Te0_err = np.asarray(c.get(r'\results::thomson.profiles.auto:te:error_bar'));
rhoT = np.asarray(c.get(r'\results::thomson:profiles:auto:rho')); #     TCV standard poloidal rho=sqrt(psi_bar)= linspace(0,1,41)
tT = np.asarray(c.get(r'\results::thomson.profiles.auto:time'));# s   time vector (usually shorter than time vector of raw data) ne.dim{1};

time = 0.9
it2 = np.argmin(abs(tvec-time))
time = tvec[it2]
it1 = np.argmin(abs(tT -time))

from matplotlib.pylab import *

#import IPython
#IPython.embed()

def plot_thompson(t):
    it2 = np.argmin(abs(tvec-t))
    time = tvec[it2]
    it1 = np.argmin(abs(tT -t))
    subplot(121)
    title(t)

    errorbar(psi[:,it2],ne[:,it2],ne_err[:,it2],fmt='o')
    fill_between(rhoT,ne0[:,it1]-ne0_err[:,it1],ne0[:,it1]+ne0_err[:,it1],
                alpha=.2, facecolor='b',edgecolor='None')
    plot(rhoT,ne0[:,it1])
    xlim(0,1)
    ylim(0,ne.max())
    subplot(122)
    title(t)

    errorbar(psi[:,it2],Te[:,it2],Te_err[:,it2],fmt='o')
    fill_between(rhoT,Te0[:,it1]-Te0_err[:,it1],Te0[:,it1]+Te0_err[:,it1],
                alpha=.2, facecolor='b',edgecolor='None')
    plot(rhoT,Te0[:,it1])
    xlim(0,1)
    ylim(0,Te.max())
    draw()
    
plot_thompson(0.6);
show()

    

fig = figure(figsize=(10,5))
fig.show()
for t in  tT:
    fig.clf()

    pause(.1)
    
    
    
    
    
    



