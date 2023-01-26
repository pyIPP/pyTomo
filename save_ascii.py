import numpy as np
import sys,os

tmp_path = os.path.expanduser('~/tomography/tmp/')
files = os.listdir(tmp_path)

for f in files:
    if f.startswith('postprocessing'):
        data = np.load(tmp_path+f, allow_pickle=True)
        r = data['rho']
        t = data['tvec']
        emiss = data['fsa_emiss']
        header = 'time [s]\n'
        header += ','.join([str(tt) for tt in np.round(t,4)])
        header += '\nrho\n'
        header += ','.join([str(rr) for rr in np.round(r,4)])
        header += '\nTotal radiation [Wm^-3]\n'
        print('saving to '+ f[:-4]+'.txt')
        np.savetxt(f[:-4]+'.txt',emiss,  header = header, fmt='%5.3e')



    if f.startswith('Emissivity_'):
        data = np.load(tmp_path+f, allow_pickle=True)
        t = data['tvec']
        r = data['rvec']
        z = data['zvec']
        g = data['gres']
        g = g*data['gres_norm']
        
        header = 'time [s]\n'
        header += ','.join([str(tt) for tt in np.round(t,4)])
        header += '\nR [m]\n'
        header += ','.join([str(rr) for rr in np.round(r,4)])
        header += '\nZ [m]\n'
        header += ','.join([str(zz) for zz in np.round(z,4)])
        header += '\nradiation [Wm^-3] (data are flattened in MATLAB order)\n'
        np.savetxt(f[:-4]+'.txt',g.flatten(order='f'),  header = header, fmt='%5.3e')


