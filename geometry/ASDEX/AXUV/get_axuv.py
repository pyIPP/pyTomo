    
    
#from loader import * 
#from scipy.interpolate import interp1d
#import os,sys
from collections import OrderedDict

#from matplotlib.pylab import *
import time 
from numpy import *

    
    

    
class loader_AXUV():
    
    
    """
    class for the loading of the AXUV data from AUG

    """


    n_smooth = 16

    
    def __init__(self, shot, geometry_path,dd,experiment='AUGD',edition=0):
   
        """

        :var int shot:  Number of selected shot
        :var str geometry_path: Path to saved geometry, boundary or data
        :var object dd: loading library for AUG data
        :var str experiment: name of AUG  experiment
        :var int edition: edition of AUG  experiment

        """
      
        
        self.geometry_path = geometry_path

        self.shot = shot
        self.experiment = experiment
        self.ed = edition
        self.dd = dd
        self.wrong_dets_damaged = []#BUG 

       
        calib_shot = self.dd.cShotNr('AUGD', 'BLC',self.shot)
       
        if not self.dd.Open( 'BLC', calib_shot):
            #second try
            calib_shot = self.dd.cShotNr('AUGD', 'BLC', calib_shot-1)
            if not self.dd.Open( 'BLC', calib_shot):
                raise Exception('No BLC shotfile!')
            
       


        names = [n.strip() for n in self.dd.GetNames()]
        self.all_names = [n for n in names if n[:2] == 'CS']

                    

        self.calib  = {}
        self.offset = {}

        self.channels= {}
        self.Cams = {}
        self.active = {}

        for name in self.all_names:
            ncalib = int(self.dd.GetParameter(name, 'NCALSTEP' ))
            multi = []
            offset = 0
            for i in range(ncalib):
                multi.append(self.dd.GetParameter(name, 'MULTIA%.2d'%i))
                offset = offset*multi[-1]+self.dd.GetParameter(name,'SHIFTB%.2d'%i)
            self.calib[name[1:]]  = prod(multi)
            self.offset[name[1:]] = offset
            channel = self.dd.GetParameter(name, 'Channel')
            camera = self.dd.GetParameter(name, 'Cam').item().strip().decode('utf-8')

            if channel!= 0:
                if not camera in self.Cams:
                   self.Cams[camera] = []
                self.Cams[camera].append(name[1:])
                self.channels[name[1:]] = channel




        all_cams = unique(list(self.Cams.keys()))
        active_los = {s:bool_(self.dd.GetParameter(s, 'active' )) for s in all_cams}
        activ_cam = [c for c in all_cams if any(active_los[c]) and not c in ('DT1', 'DT3')]
        nlos =  [sum(active_los[c]) for c in activ_cam]
        ind = argsort(nlos)
        activ_cam = [activ_cam[i] for i in ind[::-1]]
      


        self.Phi     = {s:self.dd.GetParameter(s, 'P_Blende' )[active_los[s]] for s in activ_cam}
        self.R_start = {s:self.dd.GetParameter(s, 'R_Blende' )[active_los[s]] for s in activ_cam}
        self.z_start = {s:self.dd.GetParameter(s, 'z_Blende' )[active_los[s]] for s in activ_cam}
        self.R_end   = {s:self.dd.GetParameter(s, 'R_end'    )[active_los[s]] for s in activ_cam}
        self.z_end   = {s:self.dd.GetParameter(s, 'z_end'    )[active_los[s]] for s in activ_cam}
        self.theta   = {s:arctan2(self.z_end[s]-self.z_start[s],self.R_end[s]-self.R_start[s])  for s in activ_cam}
        self.delta   = {s:self.dd.GetParameter(s, 'delta'    )[active_los[s]] for s in activ_cam}
        self.sig_dict= {s:[c.decode('utf-8') for c in self.dd.GetParameter(s,'RAW')] for s in activ_cam}
        self.dd.Close()
        #import IPython
        #IPython.embed()
  
        
        self.subcam_ind = {c:[self.delta[c]+self.R_start[c] == r for r in unique(self.delta[c]+self.R_start[c])] for c in activ_cam}
        self.DAS_limit = 2**14-1


        #geometric corrections  only estimated!!!
        n =  self.R_end['DHC']-self.R_start['DHC'], self.z_end['DHC']-self.z_start['DHC']
        corr = 0.04#rad
        self.R_end['DHC'][self.subcam_ind['DHC'][1]]+= -n[1][self.subcam_ind['DHC'][1]]*tan(corr)
        self.z_end['DHC'][self.subcam_ind['DHC'][1]]+= +n[0][self.subcam_ind['DHC'][1]]*tan(corr)
            
        n =  self.R_end['DVC']-self.R_start['DVC'], self.z_end['DVC']-self.z_start['DVC']
        corr = 0.02#rad
        self.R_end['DVC'][self.subcam_ind['DVC'][0]]+= -n[1][self.subcam_ind['DVC'][0]]*tan(corr)
        self.z_end['DVC'][self.subcam_ind['DVC'][0]]+= +n[0][self.subcam_ind['DVC'][0]]*tan(corr)
        
        self.activ_cam = activ_cam

         
        self.scam_ind = OrderedDict()
        self.cam_ind = OrderedDict()
        self.detectors_dict = OrderedDict()
        
          
        nlos = 0
        self.all_los = []
        for c in activ_cam:
            for isub,ind in enumerate(self.subcam_ind[c]):
                name = '' if len(self.subcam_ind[c])==1 else '_'+str(isub+1)
                self.scam_ind[c+name] = where(ind)[0]+nlos
            self.cam_ind[c] = nlos+arange(len(ind))
            nlos+= sum(active_los[c])
            self.detectors_dict[c] = [c+'_%d'%self.channels[n] for n in self.Cams[c]]
            self.all_los.extend(self.detectors_dict[c])


        self.nl = nlos
        self.activ_subcam = self.scam_ind.keys()

        self.calb_0 = ones(len(self.activ_subcam))
        self.dets = arange(self.nl)
        
        
        
        
        self.dets_index = [v for k,v in self.cam_ind.items()]

        #separate for each subcamera
        ##self.dets_index = [v for k,v in self.scam_ind.items()]
        ##self.detectors_dict = self.scam_ind


     
        self.dd.Open('XVR', self.shot)
        tlen = self.dd.GetInfo('Dio-Time').tlen
        tbeg = squeeze(self.dd.GetTimebase('Dio-Time',1,1))
        tend = squeeze(self.dd.GetTimebase('Dio-Time',tlen,tlen))
        self.tvec = linspace(tbeg, tend,tlen) #BUG assume equally spaced time vector, but faster loading, double accuracy
        self.dd.Close()
        
        

    
    def get_data(self,tmin=-infty,tmax=infty):
        """ Main diagnostics are diod Bolometers,
        """

        bolo_shotfiles = 'XVR', 'XVS'

        t_offset = -0.1
        nbeg,nend,noffset = self.tvec.searchsorted([tmin,tmax,t_offset])
        nend = min(nend+self.n_smooth//2, len(self.tvec)-1)
        nbeg = max(nbeg-self.n_smooth//2,0)
        #problems with numerical precision of the singles
        tvec = self.tvec[nbeg:nend]
        nt = nend-nbeg
       
        
        #load raw data
        data = empty((nt,self.nl),dtype=int16)
        offset = zeros((noffset, self.nl),dtype=int16)
        calib = zeros(self.nl)
        for bs in bolo_shotfiles:
            #if bs == bolo_shotfiles[1]:
            self.dd.Open(bs, self.shot)
                
            names = [n.strip() for n in self.dd.GetNames() if n[0]=='S']
            for cam in self.activ_cam:
                ind = self.cam_ind[cam]
                signals = self.sig_dict[cam]
                for ilos,sig in zip(ind,signals):
                    sig = sig.strip()
                    if sig in names:
                        data[:,ilos] = self.dd.GetSignal(sig,False, nbeg+1,nend) 
                        offset[:,ilos] = self.dd.GetSignal(sig,False,1,noffset)
                        calib[ilos] = self.calib[sig]
            self.dd.Close()

      
        #averadge over "n_smooth" points 
        data = mean(data[:(nt//self.n_smooth)*self.n_smooth].reshape(nt//self.n_smooth, self.n_smooth,self.nl),1,dtype='float32')
        tvec = mean(tvec[:(nt//self.n_smooth)*self.n_smooth].reshape(nt//self.n_smooth, self.n_smooth),1)

        mad = lambda x: mean(abs(x-mean(x,0)),0)*1.4
        noise = mad(offset)
        offset = mean(offset,0)
        data -= offset
        data *= calib
        das_upper_limit = calib*(self.DAS_limit-offset)
        das_lower_limit = -calib*(offset)
        noise*= calib

        data_err = ones_like(data)
        data_err*= noise
        data_err[(data>das_upper_limit*0.99 )|(data<=das_lower_limit*0.99)] = infty  
        
    
        
        if self.shot > 31591:
            pass
        elif self.shot > 30135:
            pass
        elif self.shot > 28523:
            pass      
        elif self.shot > 27352:
    
            ind = self.scam_ind['DHC_2']
            data[:,ind]/= 1.05
            
            ind = self.scam_ind['DHC_3']
            data[:,ind]*= 1.3
            
            ind = self.scam_ind['DHC_1']
            data[:,ind]*= 1.1


        return tvec, data, data_err

 
 
    def load_geom(self, path):

        for i,det in enumerate(self.activ_cam):
            R1 = self.R_start[det]
            Z1 = self.z_start[det]
            R2 = self.R_end[det]  
            Z2 = self.z_end[det]  
            savetxt(self.geometry_path+'/detector_%s_x.txt'%det,c_[R2,R1],fmt='%.4f')
            savetxt(self.geometry_path+'/detector_%s_y.txt'%det,c_[Z2,Z1],fmt='%.4f')

     
     
 

def main():
    import os
    import os,sys
    sys.path.append('/afs/ipp/home/t/todstrci/TRANSP/')
    
    axuv = loader_AXUV(30382)
    T = time.time()
    axuv.get_data(4.6, 4.7)
    
    print( 'loaded in ',time.time()-T)
    axuv.load_geom()

    #tvec, data, error = sxr.get_signal_fast(4.6, 4.7)
    
    #sxr.hardcoded_corrections(tvec, data)
    #print 'finished'
 
    #ind = argmax(data.max(0))
    #errorbar(tvec, data[:,ind], error[:,ind])
    #show()
    
    
    
    
    

if __name__ == "__main__":
    main()
