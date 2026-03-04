#!/usr/bin/python
# -*- coding: utf-8 -*-
# ##############################################################################################
# ##############################################################################################
#
#  Arguably high-level class to access JET PPF data
#  
#  24/09/10 : Initial version, TG
#
# Contact : thomas.gerbaud@jet.uk
#
# ##############################################################################################
# ##############################################################################################
 
import warnings
warnings.filterwarnings('ignore', '.*')
 
# from include import  *
import os
import sys
import time
from copy import deepcopy as dcopy
 
import numpy as N
 
from numpy import *
 
from scipy import vstack, hstack, tile, interpolate
 
from multiprocessing import Pool, TimeoutError, Process, Lock
 
 
# separate utilities file
from utils import unique,interp,whoisdad,whoami,set_ax,smooth
 
 
 
# ##############################################################################################
# ##############################################################################################
#  Available PPF signals 
#  -  "dda" is the dda
#  -  "data" are the dtyp you want to read. Must usually have the same time dimensions
#  -  "signals" the name of each data you read - I feel it more convenient 
#  -  "units" is the actual PPF unit for each data
#  -  "data", "signals" and "units" are ordered lists 
#  -  time vectors are automatically created and added to "signals"
# ##############################################################################################
# ##############################################################################################
 
USE_THREADING = True
N_PROC = 10
SIG_THREADING = ['KK3 - full']
 
SIG_LIST = { 
 
    # #####################################################################################
    # KK1 - Michelson Interferometer Te 
    # #####################################################################################
 
    'KK1':{'dda':'ECM1', 'data':['PRFL'], 'signals':['Te', 'r'], 'units':['eV']}, 
 
    'KK1 - Tmax':{'dda':'ECM1', 'data':['TMAX'], 'signals':['Tmax'], 'units':['eV']}, 
 
    # #####################################################################################
    # KK3 - ECE radiometer - Te
    # #####################################################################################
 
    # low time dynamic
    'KK3 - short':{'dda':'KK3', 'data':['TPRF', 'CPRF'], 'signals':['Te', 'r'], 'units':['eV']}, 
 
    # high time dynamic
    'KK3 - full': {'dda':'KK3', 
                   'data':[['TE%.2i'%x for x in arange(1,97)], ['RC%.2i'%x for x in arange(1,97)]],
                   'signals':['Te', 'r'], 'units':['keV']},
 
    # #####################################################################################
    # KE3 - CORE LIDAR - Te on magnetic axis 
    # #####################################################################################
 
    # Te on magnetic axis 
    'KE3 - Te0':{'dda':'LIDX', 'data':['TE0'], 'signals':['Te'], 'units':['eV']}, 
 
    # Te
    'KE3':      {'dda':'LIDR', 'data':['TE','NE'], 'signals':['Te', 'Ne', 'r'], 'units':['eV','m^-3']}, 
 
 
    # #####################################################################################
    # KE9 - EDGE LIDAR
    # #####################################################################################
 
    # not good R position ...
    'KE9':{'dda':'KE9D', 'data':['TE'], 'signals':['Te', 'r'], 'units':['eV']}, 
 
 
    # #####################################################################################
    # KE11 - HRTS - Te, Ne
    # #####################################################################################
 
    'KE11':{'dda':'HRTS', 'data':['TE', 'NE'], 'signals':['Te', 'Ne', 'r'], 'units':['eV','m^-3']}, 
 
    # #####################################################################################
    # KG1V - interferometry - Ne
    # #####################################################################################
 
    'KG1V':{'dda':'KG1V', 'data':['LID3', 'LID4'], 'signals':['Ne - LID3', 'Ne - LID4'], 'units':['m^-2', 'm^-2']}, 
 
    # #####################################################################################
    # KG10 - reflectometry - Ne
    # #####################################################################################
 
    'KG10':{'dda':'KG10', 'data':['NE', 'R'], 'signals':['Ne', 'r'], 'uid':'chain1', 'units':['m^-3']}, 
 
    # #####################################################################################
    # MAGN - Ip, Bvac
    # #####################################################################################
 
    'MAGN':{'dda':'MAGN', 'data':['IPLA', 'BVAC'], 'signals':['Ip', 'BVAC'], 'units':['-A', '-T']},
 
    # #####################################################################################
    # Additional power
    # #####################################################################################
 
    'ICRH':{'dda':'ICRH', 'data':['PTOT'], 'signals':['PTOT']},
    'LHCD':{'dda':'LHCD', 'data':['PTOT'], 'signals':['PTOT']},
    'NBI': {'dda':'NBI',  'data':['PTOT'], 'signals':['PTOT']},
 
    # #####################################################################################
    # EFIT - Equilibrium
    # #####################################################################################
 
    # B field
    'EFIT_B': {'dda':'EFIT', 'data':['BRAX','BTAX','BZAX' ], 'signals':['BR', 'BT', 'BZ', 'r']},
 
    # Boundary thingsx
    'EFIT_FBND':   {'dda':'EFIT', 'data':[ 'FBND'], 'signals':[ 'FBND']},
    'EFIT_RBND': {'dda':'EFIT', 'data':[ 'RBND'], 'signals':[ 'RBND', 'r']},
    'EFIT_ZBND': {'dda':'EFIT', 'data':[ 'ZBND'], 'signals':[ 'ZBND', 'z']},
 
    # #####################################################################################
    # Dalpha 
    # #####################################################################################
 
    'Da': {'dda':'S3AD', 'data':[ 'AD35'], 'signals':[ 'Da']},
    
   
    
    'BOLO': {'dda':'BOLO', 'data':['ETEN'], 'signals':['ETEN']}
 
}
 
 
 
 
# ##############################################################################################
# ##############################################################################################
# Useful programs 
# ##############################################################################################
# ##############################################################################################
 
 
# ##############################################################
# read a value with ppfget
# ##############################################################
 
def getppf(pulse, dda, dtyp, uid = 'JETPPF', seq = 0,  verbose = False):
 
    sys.path.append('/jet/share/lib/python')
    from ppf import ppfgo, ppfget, ppfuid, ppfclo
 
    default = [False]*3
    def close_and_return(x):
        ppfclo(pulse, dda, seq)
        return x
 
    # init
    ier = ppfgo(pulse, seq = seq)
 
    if ier not in [0, 100000]:
        if verbose: 
            print(("ERROR - getppf - %i - error code %i"%(pulse, ier)))
        return close_and_return(default)
 
    # set uid
    ppfuid(uid,'r')
 
    # get data
    ihdat,iwdat,data,x,t,ier = ppfget(pulse, dda, dtyp)
 
    if ier != 0 :
        if verbose: print(('error #%i ! from (%s)'%(ier, whoisdad())))
        return close_and_return(default)
 
    if verbose:  
        print(('[PPF] : %s/%s OK'%(dda, dtyp)))
        sys.stdout.flush()
 
    # reshape
    if len(x) > 1:
        data,x,t =  reshape(data, (len(t), len(x))).T, array(x), array(t)
    else:
        data,x,t =  list(map( array, [data, x, t]))
 
    return close_and_return([data, x, t])
 
 
 
# ##############################################################################################
# ##############################################################################################
# Main class 
#   - data are read once, then stored for further use
# ##############################################################################################
# ##############################################################################################
 
class JET:
    def __init__(self, pulse=False, uid = 'jetppf', verbose = True):
        if not pulse:
            pulse = eval(input('Shot number ? : '))
 
        self.pulse = pulse
        self.uid = uid
        self.verbose = verbose
        self.cache = {}
        self.signals_with_fixed_units = []
 
    # ##############################################################################################
    #
    # ##############################################################################################
 
    def print_me(self, x):
        if self.verbose: print(x)
 
 
    # ##############################################################################################
    #
    # ##############################################################################################
 
    def print_error(self, x):
        print(('[%i] ERROR %s in %s'%(x, whoisdad())))
 
    # ##############################################################################################
    # Overwriting the [...] function.
    # ##############################################################################################
 
    def __getitem__(self, args):
        if type(args) in [list, tuple]:
            what, sig = args[0], args[1]
        else:
            what, sig = args, False
 
        if sig: 
            return self.get_signal(args[0], args[1])
 
        sig = self.main_sig(what)
        if sig is not False:
            return self.get_signal(args, sig)
 
        return self.get_signal(args)
 
 
    # ##############################################################################################
    # Overwriting the self[x]=val function.
    # ##############################################################################################
 
    def __setitem__(self, key, item): 
        if type(key) not in [list, tuple] or len(key) != 2: 
            self.print_error('need 2 arguments')
 
        what, sig = key
 
        if self.in_cache(what) and sig in self.cache[what]:
            self.cache[what][sig] = item
 
 
    # ##############################################################################################
    # Is key in cache ?
    # ##############################################################################################
 
    def in_cache(self, sig):
        return sig in self.cache
 
 
    # ##############################################################################################
    # 
    # ##############################################################################################
 
    def store_in_cache(self, sig, key, data):
        if not self.in_cache(sig):
            self.cache[sig] = {}
            for x in SIG_LIST[sig]['signals']+['t']: 
                self.cache[sig][x] = False
        self.cache[sig][key] = data
 
 
    # ##############################################################################################
    # Simple wrapper to the out-of-class getppf function  --  needed when THREADING
    # ##############################################################################################
 
    def from_ppf( self, args ):
        pulse, dda, x, uid, verbose = args
        return x, getppf( pulse, dda, x, uid=uid, verbose=verbose)
 
 
    # ##############################################################################################
    # Read parameters from SIG_LIST
    # ##############################################################################################
 
    def get_info( self, what ):
 
        sigs = SIG_LIST[what]
        dda, dtypes, signals = [dcopy(sigs[x]) for x in ['dda', 'data', 'signals']]
 
        units, uid = False, self.uid
        if 'units' in sigs: units = dcopy(sigs['units'])
        if 'uid' in sigs:   uid = dcopy(sigs['uid'])
 
        return dda, dtypes, signals, units, uid
 
    # ##############################################################################################
    # If only one data signal is read, returns it
    # ##############################################################################################
 
    def main_sig( self, what ):
        sigs = self.get_info(what)[2]
 
        if sigs.count('r'): sigs.remove('r')
        if sigs.count('t'): sigs.remove('t')
 
        if len(sigs) == 1: return sigs[0]
        return False
 
    # ##############################################################################################
    # Main function : return a signal DDA/DTYP with what=DDA and mysig=DTYP
    # If not DTYP provided it returns all the data in SIG_LIST for that DDA
    # Try to read from cache, then get hit the PPF system
    #
    # Comments : 
    #   - lots of lines are actually only used to process KK3 data.
    #   - will use multiple simultaneaous access (THREADING) for KK3 PPF as 192 DTYP are needed
    #   - that function can be include in a THREAD'D process if a Lock is provided
    # ##############################################################################################
 
    def get_signal(self,  what, mysig = False, lock=False ):
 
        def protect(f, *args):
            if lock : lock.acquire()
            out = f(*args)
            if lock : lock.release()
            return out
 
        if not self.in_cache(what):
 
            # dda, dtypes, signals, units, uid = self.get_info(what)
            dda, dtypes, signals, units, uid = protect( self.get_info, what)
 
            data, t = {}, {}
 
            # one or more signals ?
            if type(dtypes[0]) == str : mode = 0
            if type(dtypes[0]) == list: mode = 1
 
            # ##################################################################
            # get all data
            # ##################################################################
 
            # threading ? 
            use_threading =  USE_THREADING and what in SIG_THREADING
            if use_threading: 
                pool = Pool(processes=N_PROC)
                TIMEOUT = 30 #s
 
            for sig, dtyp in zip( signals, dtypes):
                if mode == 0 : dtyp_list = [dtyp]
                if mode == 1 : dtyp_list =  dtyp
 
                # ##################################################################
                if not use_threading:
                    for x in dtyp_list:
                        data[x], tmp, t[x] = self.from_ppf((self.pulse, dda, x, uid, self.verbose))[1]
 
                    # common situation : 'r' (or signals[1]) is tmp !
                    if mode == 0 and len( dtypes ) == len(signals) - 1 and signals[len(dtypes)] == 'r':
                        dtypes += [dtypes[0]+'_R']
                        data[dtypes[-1]] = tmp
 
                # ##################################################################
                if use_threading: 
 
                    args = [( self.pulse, dda, x, uid, self.verbose) for x in dtyp_list]
                    it = pool.imap_unordered( self.from_ppf, args)
 
                    while True:
                        try:
                            out  = it.next( TIMEOUT )
                        except StopIteration:
                            break
                        except TimeoutError:
                            print('Skipping (timeout)')        
                            continue
 
                        x = out[0]
                        data[x], tmp, t[x] = out[1] 
 
            # ##################################################################
            # \ get all data
            # ##################################################################
 
            # do we have some data ?
            if all( [x is False for x in list(data.values())]):
                self.cache[what] = False 
                return False
 
            # check time - all signals must have the same shape
            if mode == 0:
                for x in [data, t]:
                    if len(x) > 1 and len(unique(list(map(len, list(x.values()))))) > 1:
                        self.print_error('while reading %s : dimensions mismatch'%(dda))
                        return False
 
            # ( usually skipped ) 
            # check time - all signals from a same dtype must have the same shape
            if mode == 1:
                for x in [data, t]:
                    for dtyp in dtypes:
                        tmp = [len(x[y]) for y in dtyp]
                        if len(tmp) > 1 and len(unique(tmp)) > 1:
                            self.print_error('while reading %s : dimensions mismatch for %s...%s'%(dda, dtyp[0], dtyp[-1]))
 
            # ( usually skipped ) 
            # reshape data if needed   
            if mode == 1:
                # reshape
                for sig, dtyp in zip(signals, dtypes):
                    t[sig] = t[dtyp[0]] 
                    data[sig] = zeros((len(dtyp), len(data[dtyp[0]])))
                    for ii,dd in enumerate(dtyp):
                        data[sig][ii,:] = data[dd]
                        data.pop(dd)
                        t.pop(dd)
 
 
            # ( usually skipped ) 
            # interp data if needed   
            if len(unique(list(map(len, list(t.values()))))) > 1:
                sig_ref = signals[argmax([len(t[x]) for x in signals])]
                t_ref = t[sig_ref]
                for sig in signals:
                    if sig != sig_ref:
                        ind = [jj for jj,val in enumerate(t_ref) if val>=t[sig][0] and val<=t[sig][-1]]
                        x = zeros((data[sig].shape[0], len(t_ref)))
                        for ii,dat in enumerate( data[sig]):
                            x[ii,ind] = interp( t[sig], dat, t_ref[ind], 'linear')
 
                        t[sig], data[sig] = x, t_ref
 
            # save in cache
            if mode == 0: keys = dtypes
            if mode == 1: keys = signals
 
            for sig,key in zip(signals,keys):
                 # self.store_in_cache(what, sig, data[key])  
                protect( self.store_in_cache, what, sig, data[key])  
 
            # self.store_in_cache(what, 't', t.values()[0])
            protect( self.store_in_cache, what, 't', list(t.values())[0])
            # self.fix_units( what )
            protect( self.fix_units, what )
        #
 
        if mysig is False:
            return self.cache[what]
        else:
            if self.cache[what] is not False:
                return self.cache[what][mysig]
            else:
                return False
 
 
 
    # ##############################################################################################
    # *args version
    # ##############################################################################################
    def get_signal_with_args(self,  args):
        return self.get_signal(*args)
 
    # ##############################################################################################
    # Change units (once):
    #    - eV              ->             keV
    #    - m^-2            ->         1e-20*m^-2 
    #    - m^-3            ->         1e-19*m^-3  
    #    - A               ->    (opposite sign) kA
    #    - T               ->    (opposite sign) T
    # ##############################################################################################
 
    def fix_units(self, what):
        dda, dtypes, signals, units, uid = self.get_info(what)
 
        if units is False: return
 
        # already fixed ?
        fixed = self.signals_with_fixed_units
        if fixed.count( what ): 
            return None
        else:
            fixed.append(what)
 
        for sig in ['r', 't']:
            if signals.count(sig): 
                signals.remove(sig)
 
        def f(sig, myunit, from_unit, to_unit, corr ):
            if myunit == from_unit:
                self.cache[what][sig] *= corr
                return to_unit
                # self.print_me('-> [%s,%s] modified to %s by x%f'%(sig, myunit, unit, corr))
 
 
        for ii in arange(len(signals)):
            sig, unit = signals[ii], units[ii]
 
            units[ii] = f( sig, unit,   'eV',        'keV',   1e-3)  
            units[ii] = f( sig, unit, 'm^-3',  '1e19 m^-3',  1e-19)  
            units[ii] = f( sig, unit, 'm^-2',  '1e20 m^-2',  1e-20)  
            units[ii] = f( sig, unit,   '-A',         'kA',  -1e-6)  
            units[ii] = f( sig, unit,   '-T',          'T',     -1)     
 
        SIG_LIST[what]['units'] = units
 
 
    # ##############################################################################################
    # find Tmax for the pulse
    # ##############################################################################################
 
    def T_max(self):
        if not self.in_cache('T max'):
            T1, T2 = self['KK1 - Tmax'], self['KE3 - Te0'] 
            T1, T2 = [smooth(x, 10) for x in [T1, T2]]
 
            self.cache['T max'] = .5*(T1.max()+T2.max())
        #
 
        return self.cache['T max']
 
    # ##############################################################################################
    # find pulse start/end (based on Ip)
    # ##############################################################################################
 
    def t_pulse(self):
        if not self.in_cache('t pulse'):
            Ip, t = self['MAGN', 'Ip'], self['MAGN', 't']
            ind = where( Ip > 0.05)[0] # 50kA
 
            self.cache['t pulse'] = [t[ind][0], t[ind][-1]]
        #
 
        return self.cache['t pulse']
 
 
    # ##############################################################################################
    # Easy plot
    # ##############################################################################################
 
    def plot_signal(self, what, sig=False, ifig=False):
 
        r = False
        data = self[what]
 
        if data is False:
            return 
 
        if sig is False:
            sig = self.main_sig(what)
 
        if type(data) != dict: 
            t = self[what, 't']
            try:
                r = self[what, 'r']
            except:
                pass
        else:
            if 'r' in data:
                data, r, t = [data[x] for x in [sig, 'r', 't']]
            else:
                data, t = [data[x] for x in [sig, 't']]
 
        from pylab import figure
 
        if ifig is False: ifig = 0 
        f = figure(0)
        f.clf()
        ax = f.add_subplot(111)
 
        if what in ['KG10']:
            data, r, t = data[:, ::100], r[:, ::100], t[:,::100]
 
        sigs = SIG_LIST[what]
        units = sigs['units'][ sigs['signals'].index(sig)]
 
        t_ok, T = self.t_pulse(), self.T_max()
        if r is not False:
            ax.plot( r, data, 'k.', ms=10, alpha=.5)
            xl = 'r [m]'
            ax.axis([1.8, 4.0, 0, T])
        else:
            ax.plot( t, data, 'k-')
            xl = 't [s]'
            ax.axis((t_ok[0], t_ok[-1], 0, T))
 
 
        set_ax(ax, t='#%i - %s - %s'%(self.pulse, what, sig), xl=xl, yl=units, t_fs=16, xl_fs=14, yl_fs=14) 
        f.subplots_adjust(left = 0.1, right=0.95, bottom = 0.1, top = 0.9, wspace=0.2, hspace=0.25)
 
        f.canvas.draw_idle()
 
 
    # ##############################################################################################
    # Load common signals - try to use THREADING for quicker access
    # ##############################################################################################
 
    def load_common_signals(self, siglist = [],  use_threading=True):
 
	if siglist == [] :
	    #siglist = [   'KK1', 'KK1 - Tmax', 'KK3 - short', # Te - ECE
			#'KE3', 'KE3 - Te0', 'KE11',         # Te, ne - TS
			#'KG1V', 'KG10',                     # ne
			#'ICRH', 'NBI', 'LHCD',              # POWER
			#'MAGN', 'EFIT_B', 'Da',
			#'SXRT',
			#]
	    siglist = [ 'BOLO' ]
	
        use_threading = use_threading  and  USE_THREADING 
 
        # ##################################################################
        if not use_threading:
            list(map( lambda x:self[x], siglist))
 
        # ##################################################################
        if use_threading: 
            lock = Lock()
 
            plist = []
            for sig in siglist:
                plist.append( Process(target=self.get_signal, args=(sig, False, lock,)))
                plist[-1].start()
 
            for p in plist:
                p.join( 30 ) # 30s timeout
 
        # finish
        self.print_me('All data read')
   
  
def datasend(datafile):
    try:
	print('sending ...')
	from ftplib import FTP
	ftp = FTP('193.86.154.26')
	psw = 'YWJjKzEyMw==\n'
	ftp.login('JET', psw.decode('base64'))
	ftp.cwd('/tmp/harddisk/JET')
	#ftp.retrlines('LIST')
	ftp.storbinary('STOR '+datafile, file(datafile, 'rb'))
	ftp.quit()
	print('data send', datafile)
    except:
	raise
	print('connection lost', datafile)
 

from numpy import *
from scipy.io import loadmat, savemat
import multiprocessing

def main():
    UPLOAD=False
    pulselist =  [65670]
    for PULSE  in pulselist:
	siglist = [ 'BOLO' ]
	eq = JET(PULSE , verbose=True)
	eq.load_common_signals(siglist)
	for SIG in siglist:
	    datafile = 'DATA_'+str(PULSE)+'_'+str(SIG)
	    save(datafile, eq[SIG])
	    savemat(datafile, {'data': eq[SIG]})

	    if UPLOAD:
                #p = multiprocessing.Process(target=datasend, args=(geometry_path, 'Data_'+str(pulse)+'_'+SOURCE+'.npz'))
                #p = multiprocessing.Process(target=datasend, args=(datafile+'.npy'))
                #p.start()
                #sleep(10)
                #p2 = multiprocessing.Process(target=datasend, args=(datafile+'.mat'))
                #p2.start()
                datasend(datafile+'.npy')
                datasend(datafile+'.mat')
	   

    #T,r,t = eq['KK1'], eq['KK1', 'r'], eq['KK1', 't']  # Retrieve Te(r,t) from KK1
    #SXRT = eq['SXRT']
    #print SXRT
    #save('SXRT', SXRT)
    
    
if __name__ == "__main__":
    main()
    