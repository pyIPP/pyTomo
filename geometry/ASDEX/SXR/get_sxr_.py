
#from loader import * 
from scipy.interpolate import interp1d
import os,sys
from collections import OrderedDict
from matplotlib.pylab import *
import time 
from numpy import *
#import numpy as np

from shared_modules import fast_svd





class loader_SXR():
    
    
    """
    class for the loading of the SXR  data from AUG

    """
    

    sigma = 0.02
    def __init__(self, shot,geometry_path,dd,fast_data=False,  experiment='AUGD',edition=0):
       
        """

        :var int shot:  Number of selected shot
        :var str geometry_path: Path to saved geometry, boundary or data
        :var object dd: loading library for AUG data
        :var str experiment: name of AUG  experiment
        :var int edition: edition of AUG  experiment

        """
        
        self.shot = shot
        self.experiment = experiment
        self.fast_data = fast_data
        self.ed = edition
        self.geometry_path = geometry_path
        self.dd = dd
        self.calib_shot = self.dd.cShotNr('AUGD', 'CSX',self.shot)
        self.dd.Open( 'CSX', self.calib_shot)
        names = self.dd.GetNames()  

        
        signals = [b[1:] for b in names if b[0]== 'C']
 
        
        shotfiles      = [  self.dd.GetParameter('C'+s, 'SX_DIAG'  ).tostring() for s in signals]
        geometry_names = ('Tor_Pos','RPINHOLE','ZPINHOLE','REND','ZEND','D_Width',
                    'D_Length','THETA','CAMANGLE','D_Gap','P_Width','P_Length','Foc_Len')
                    
        self.geometry = {}
        for name in geometry_names:
            self.geometry[name] = {s:self.dd.GetParameter('C'+s,name) for s in signals}
        
        self.status    = {s:self.dd.GetParameter('C'+s, 'ADDRESS')!=256 for s in signals}
        self.SampFreq  = {s:self.dd.GetParameter('C'+s, 'SampFreq') for s in signals}
        filt_mat       = {s:self.dd.GetParameter('C'+s, 'FILT-MAT').item() for s in signals}
        thickness      = {s:self.dd.GetParameter('C'+s, 'THICKNES') for s in signals}
        self.ADCrange  = {s:self.dd.GetParameter('C'+s, 'ADCrange') for s in signals}
        self.ADCmin = 0

        
        self.different_det = {s:(abs(thickness[s]-75e-6)>1e-5) | (filt_mat[s]!= 'Be') for s in signals}



        self.MULTIA   = {}
        self.SHIFTB   = {}

        for k,s in enumerate(signals):
            n = self.dd.GetParameter('C'+s, 'NCALSTEP')
            self.MULTIA[s] = [self.dd.GetParameter('C'+s,'MULTIA%.2d'%i) for i in range(n)]
            self.SHIFTB[s] = [self.dd.GetParameter('C'+s,'SHIFTB%.2d'%i) for i in range(n)]
    

        self.SXR_diods = OrderedDict()
        for sf in unique(shotfiles):  self.SXR_diods[sf] = []

        for sf,sig in zip(shotfiles,signals):
            if sig[0] != 'T': #remove T camera
                self.SXR_diods[sf].append(sig)
            
        self.SXR_diods.pop('OOO')  #null detector

        
        #identify a single subcameras of the SXR system
        self.all_los = hstack(self.SXR_diods.values())
        self.all_los.sort()
        
        
        #print  self.all_los[43]

        cams, index = unique([l[0] for l in self.all_los], return_inverse=True)
        CANGLE = array([self.geometry['CAMANGLE'][l] for l in self.all_los])
        self.detectors_dict=OrderedDict()
        for i,c in enumerate(cams):
            uangle,sub_ind,sub_index = unique(CANGLE[index==i], return_inverse=True, return_index=True)
            n_cams = len(uangle)
            if n_cams == 1:
                self.detectors_dict[c] = self.all_los[index == i]
            else:
                sort_ind = argsort(sub_ind)
                for j in range(n_cams):
                    self.detectors_dict[c+str(j+1)] = self.all_los[index == i][sub_index == sort_ind[j]]
        
        #find channels corresponding to the each subcamera
        self.cam_ind = OrderedDict()
        for k,item in self.detectors_dict.iteritems():
            self.cam_ind[k] = zeros(len(item), dtype=int)
            for i,los in enumerate(item):
                self.cam_ind[k][i] = where(self.all_los==los)[0][0]
  

        self.nl = len(self.all_los)
        self.calb_0 = ones(len(self.detectors_dict.keys()))
        self.Ndets = len(self.detectors_dict.keys())
        self.dets_index = [v for k,v in  self.cam_ind.iteritems()]
        self.dets = arange(self.nl)
                
       
        try:
            SXerrors = genfromtxt('/afs/ipp-garching.mpg.de/home/s/sxry/SXR/errdata/SXerrors/SX_err_%d.dat'%self.shot,
                        dtype={'names': ('channel', 'Diag', 'signal','error','13or14','15or16'),
                                'formats': ('S5', 'S3','i4','i4','i4','i4')},usecols=(0,1,2,3,4,5))
            print '\n\nSuspicious channels from /home/s/sxry/SXR/errdata/SXerrors/SX_err_%d.dat'%self.shot
            print '\t'.join(('channel', 'name','errors','13or14','15or16'))
            
            CSI="\x1B["
            reset=CSI+"m"
            red_start = CSI+"31;40m"
            red_end = CSI + "0m" 
            wrong_ch = []
            for i,los in enumerate(self.all_los):
                for err in SXerrors:
                    if los== err[0]: 
                        if err[3] < 0 or err[3] > 10000:
                            print red_start,
                            wrong_ch.append(i+1)
                        print i+1,'\t', err[0],'\t', err[3],'\t',err[4],'\t',err[5],
                        if err[3] < 0 or err[3] > 10000: print red_end
                        else: print 
                        
            print wrong_ch
         
                        
        except Exception as e:
            print e
            pass
                       
        
        #print self.all_los[[148, 149, 150,151,153, 192, 59]]
        
        #exit()
        
        #exit()
        self.Phi = [deg2rad( self.geometry['Tor_Pos'][s]+45) for s in self.all_los]#add 45deg to make tor postion consistent with diaggeom
        SampFreq = [self.SampFreq[s] for s in self.all_los]
        #savez('SXR_data', Phi=self.Phi, SampFreq=SampFreq)


        if fast_data:
            if not self.dd.Open('SXA',self.shot):
                raise Exception('SXA shotfile is not avalible')

            #trick to avoid loading the whole time base!
            tlen = self.dd.GetInfo('Time').tlen
            tbeg = self.dd.GetTimebase('Time',1,1)[0]
            tend = self.dd.GetTimebase('Time',tlen,tlen)[0]
            self.tvec = linspace(tbeg, tend, tlen)  #double precision

        else:
            if not self.dd.Open('SSX', self.shot, experiment=self.experiment, edition=self.ed):
                raise Exception('SSX shotfile is not avalible')
            self.tvec  = self.dd.GetTimebase('time' )

        self.dd.Close()
        

    def get_data(self,tmin=-infty,tmax=infty):

        if self.fast_data:
            return self.get_data_fast(tmin,tmax)


        self.dd.Open('SSX', self.shot, experiment=self.experiment, edition=self.ed)
        tvec  = self.tvec
        tvec2 = self.dd.GetTimebase('time2')

        imin,imax = tvec.searchsorted([tmin, tmax])
        if imax!= len(tvec): imax+= 1
    
        
        
        imin2,imax2 = tvec2.searchsorted([tmin,tmax])
        if imin2 == imax2: imax2+= 2
        if imax2!= len(tvec2): imax2+= 1
        

        tvec  = tvec[ imin:imax]
        tvec2 = tvec2[imin2:imax2]

        signals = self.dd.GetNames() 


        #load uncalibrated data - it is  faster!!
        data = empty((imax-imin,self.nl), dtype=int16)
        for ilos,los in enumerate(self.all_los):
            data[:,ilos] = self.dd.GetSignal(los,False,imin+1,imax)

        #print self.all_los
        #identify wrong time points 
        ADCrange = array([self.ADCrange[s] for s in self.all_los])
        wrong_poinstDAS = data > ADCrange #DAS failures 
        wrong_poinst = (data == self.ADCmin) | (data>=ADCrange)  #points out of DAS range
        
        offset = zeros(self.nl,dtype=single)
        calib  =  ones(self.nl,dtype=single)
        


        for j, name in enumerate(self.all_los):
            M = self.MULTIA[name]
            S = self.SHIFTB[name]
            offset[j] = ((S[0]*M[1]+S[1])*M[2]+S[2])*M[3]+S[3]
            calib[j] = prod(M)
        
        #calibrate signals
        data = single(data)
        data *= calib
        data += offset
        
    
        
        data_min,data_max = None,None
        if self.all_los[0]+'_MN' in signals:
            data_min = empty((imax2-imin2, self.nl),dtype=int16)
            data_max = empty((imax2-imin2, self.nl),dtype=int16)
            for ilos,los in enumerate(self.all_los):
                data_min[:,ilos] = self.dd.GetSignal(los+'_MN',False,imin2+1,imax2)
                data_max[:,ilos] = self.dd.GetSignal(los+'_MX',False,imin2+1,imax2)

        

        sigma_name = 'Sm' if not self.experiment in ["SXRY",'TODSTRCI'] and self.shot < 30619  else 'Sn'
        data_err = None
        if self.all_los[0]+'_'+sigma_name in signals:  
            data_err = empty((imax2-imin2, self.nl),dtype=int16)
            for ilos,los in enumerate(self.all_los):
                data_err[:,ilos] = self.dd.GetSignal(los+'_'+sigma_name,False,imin2+1,imax2)

        self.dd.Close()
        low_data = data
        

                
        if len(tvec)!= len(tvec2):
            trimmed_tvec = maximum(minimum(tvec[-1],tvec2),tvec[0])
            low_data = interp1d(tvec,data,axis=0,
                            kind='nearest',assume_sorted=True)(trimmed_tvec)  

        #estimate the errorbars
        if data_err is None and not data_max is None:
            #very old shotfile - naive guess of the errobars 
            data_err = abs(single(data_max-data_min))
            data_err[data_max >= ADCrange] = infty
        elif all(data_err == 0):
            #corrupted shotfile 
            data_err = data*single(self.sigma)

        elif not data_err is None and size(data_err) == self.nl*len(tvec2):

            #standard shotfile
            data_err  = single(data_err)
            data_err *= calib
            
            #compansate overetimated errobars for slow DAS
            SampFreq = array([self.SampFreq[s] for s in self.all_los])
            data_err*= SampFreq/2.e6

            if not self.experiment in ["SXRY",'MARKUSW','TODSTRCI'] and self.shot < 30619:
                #old shotfile, a different algorithm to calculate data_err was used
                data_err = minimum(data_err, (data_max-data_min)*single(calib)/2)#better error estimate in the presence of the MHD mode
                #values too close to minimum or maximum are wrong  in old shotfiles due to limited ADC range 
                data_err[((data_min == self.ADCmin)&((low_data-offset)<2*data_err))|(data_max>=ADCrange)] = infty
                            
        else:
            #when the STD is not avalible - use a guess from min and max
            data_err = single(minimum(low_data-data_min*calib-offset,data_max*calib+offset-low_data)/3)
        
        if len(tvec)!= len(tvec2):
            trimmed_tvec = maximum(minimum(tvec2[-1],tvec),tvec2[0])
            data_err = interp1d(tvec2,1./data_err/3,axis=0,bounds_error=False,
                            fill_value=infty,kind='nearest',assume_sorted=True)(trimmed_tvec)
                
    
            data_err[data_err!= 0] = 1/data_err[data_err!= 0] 
            data_err = single(data_err)
            
        data_err[wrong_poinst|isnan(data_err)|(data_err<=0)] = infty
        data[wrong_poinstDAS] = 0
        

        detector_stat = ~array([~self.status[s]|self.different_det[s] for s in self.all_los])
        useless_det = where(~detector_stat)[0]
        useless_det = self.all_los[useless_det]

        self.wrong_dets_damaged = useless_det
        self.hardcoded_corrections(tvec, data,data_err)
         
        #savez_compressed('low_tvec_'+str(self.shot), tvec=tvec, tmin=tvec[0],tmax=tvec[-1])
        #for cam, ind in self.cam_ind.iteritems():
            #print cam, data[ind].shape, data.shape
            #detector_stat = ~useless_det[ind]
            #savez_compressed('low_%s_'%cam+str(self.shot), data=data.T[ind],data_err=data_err.T[ind],data_err_disc=0*ind,detector_stat=detector_stat[ind])
        #exit()
        return tvec, data, data_err

      
        

    def hardcoded_corrections(self,tvec, data,data_err):
        
        

        
        if self.shot>=30420 and self.shot<=31900: 
            self.wrong_dets_damaged = r_[self.wrong_dets_damaged,['H_%.3d'%i for i in range(46,81)]]

        if self.shot>=30420 and self.shot<=30447: 
            self.wrong_dets_damaged = r_[self.wrong_dets_damaged,['H_%.3d'%i for i in range(0,100)]]

   


        if self.shot<=27392 and self.shot>=25995:  #wrong BE filter
            self.wrong_dets_damaged = r_[self.wrong_dets_damaged,['H_%.3d'%i for i in range(17,26)],'I_061']
        elif self.shot <= 29500:
            self.wrong_dets_damaged = r_[self.wrong_dets_damaged,['I_061',]]

            #pass
            #print '===================================================='
            print self.wrong_dets_damaged

        elif self.shot < 31776:
            self.wrong_dets_damaged = r_[['I_061',],self.wrong_dets_damaged]
            #geometry correction, 
            L_calib = r_[1.446,1.426,1.248,1.183,1.093]
            data[:,self.cam_ind['L'][:len(L_calib)]] *= L_calib

        elif self.shot < 33724:
            #etendue correction, 
            L_calib = r_[1.446,1.426,1.248,1.183,1.093]
            M_calib = r_[1.900,2.18,2.28,2.67,3.95]
            if 'L' in self.cam_ind: data[:,self.cam_ind['L'][:len(L_calib)]] *= L_calib
            if 'M' in self.cam_ind: data[:,self.cam_ind['M'][-len(M_calib):]] *= M_calib
            self.wrong_dets_damaged = r_[self.wrong_dets_damaged,['K_021',]]

        elif self.shot < 34718: 
            self.wrong_dets_damaged = r_[self.wrong_dets_damaged,
                                         ['H_048', 'H_050', 'H_052', 'H_054',
                                          'K_020', 'K_014', 'K_015', 'K_016', 
                                          'K_017', 'K_019','M_014', 'H_058','H_022']]
        elif self.shot < 34992: 
            self.wrong_dets_damaged = r_[self.wrong_dets_damaged,
                                ['K_014','K_015','K_016','K_017','K_018','K_019' ]]

        if self.shot > 33724:  #(not yet repaired)
            self.wrong_dets_damaged = r_[self.wrong_dets_damaged,['M_014','K_020','K_049','K_057','K_058']]
        
        if 34887 > self.shot > 34495: 
            self.wrong_dets_damaged = r_[self.wrong_dets_damaged,['K_053',]]
            
            

        if self.shot < 27000  and self.shot > 22000 :
            #correct the wrongly estimated offset
            if any(tvec<0.01):
                ind = any(isfinite(data_err[tvec<0.01]),axis=0)
                offset = nanmedian(data[:tvec.searchsorted(.01), ind],axis=0)
            else:
                _, data_0, data_err_0 = self.get_data(tmin=-infty,tmax=0.01)

                ind = any(isfinite(data_err_0),axis=0)
                offset = nanmedian(data_0[:, ind],axis=0)
                
            data[:,ind] -= offset

 
        
        
        
        
        


    def get_data_fast(self,tmin=-infty,tmax=infty):
        

        all_index = []
        wrong_det = []

        
        imin_,imax_ = self.tvec.searchsorted((tmin,tmax))
        tvec = self.tvec[imin_:imax_+1]

            
        ndets = 0
        isig = 0
        all_data = zeros((imax_-imin_+1, self.nl),dtype=single)
        
        
        detector_lim = zeros(self.nl)
        detector_stat = zeros(self.nl,dtype=bool)
        discr_err = zeros(self.nl)  #discretisation error
        detector_sampl = zeros(self.nl)

        

        for DAS,signals in self.SXR_diods.iteritems():
            
            
            if len(signals) == 0:
                continue

            sys.stdout.write("\rloading %2.0f%%  diag:%s  " %(isig*100./self.nl,DAS))
            sys.stdout.flush()
            
            isig+= len(signals)
            ndets+= 1
            
            
            if not self.dd.Open(DAS,self.shot):
                continue
                

            los_num = self.all_los.searchsorted(signals)
            ADCrange = self.ADCrange[signals[0]]
            SampFreq = self.SampFreq[signals[0]]

            offset = zeros(len(signals))
            calib  =  ones(len(signals))
            for j, name in enumerate(signals):
                M = self.MULTIA[name]
                S = self.SHIFTB[name]
                offset[j] = ((S[0]*M[1]+S[1])*M[2]+S[2])*M[3]+S[3]
                calib[j] = prod(M)
                
            
            #trick to avoid loading the whole time base!
            info = self.dd.GetInfo('Time')
            tlen = info.tlen
            tbeg = self.dd.GetTimebase('Time',1,1)[0]
            tend = self.dd.GetTimebase('Time',tlen,tlen)[0]
                     
            imin,imax = (r_[tmin,tmax]-tbeg)/(tend-tbeg)*(tlen-1)
            imin,imax = int(ceil(max(0,imin))), int(ceil(min(imax,tlen)))
            tbeg = self.dd.GetTimebase('Time',imin+1,imin+1)[0]
            tend = self.dd.GetTimebase('Time',imax+1,imax+1)[0]
            
            #(tmin-tbeg < 0) and (tend-tmax)>0
            nt = imax-imin+1
            data_tvec = linspace(tbeg, tend,nt )
            #import IPython
            #IPython.embed()
            
            #load data
            data = ma.zeros((nt, len(signals)),dtype=int16)
            for i,sig in enumerate(signals):
                data[:,i] = self.dd.GetSignal(sig,cal=False, nbeg=imin+1,nend=imax+1)
            self.dd.Close()
     

            ##wrong poinst, DAS failure - usually single time points 
            wrong_ind = where(data>ADCrange)
            for i,j in zip(*wrong_ind): 
                data[i,j] = 0 if i == 0 else data[i-1,j]

            
            dt = 1./SampFreq
            Nf = int(1e-3*SampFreq)

            data.mask = data <= 0  #get rid of wrong points and cut of points

            data_ = data[:(nt/Nf)*Nf].reshape(Nf,-1,size(data,1))
            data_.mask = data.mask[:(nt/Nf)*Nf].reshape(Nf,-1,size(data,1))

            #compensate the bottom cutted signal by cutting of the tops - averadge will be correct
            for i in xrange(size(data_,1)):
                for j in xrange(size(data_,2)):
                    d = data_[:,i,j]
                    if not any(d.mask) or all(d.mask):   continue
                    q = sum(d.mask)/float(size(data_,0))
                    d.mask[d >= percentile(d.data[~d.mask],(1-q)*100)] = True

     
            #overburned poinst
            data.mask[data.data == ADCrange] = False
            #negative points are corrupted measurements!
            data.mask[data.data == 0]        = False

            #apply the calibration and shift of the detectors
            data = single(data)
            data *= calib
            data += offset
            
            #interpolate on the full time resolution and fill corrupted points 
            for i,sig in enumerate(signals):
                valid = ~data.mask[:,i]
                if abs(len(data_tvec)- len(tvec)) < 0.01*len(tvec) and abs(len(data_tvec)- len(tvec))!= 0:
                    #print 
                    raise Exception('Wrong length of tvec!! %d  %d'%(len(data_tvec), len(tvec)))

                if (len(tvec) != len(data_tvec) or not all(valid)) and any(valid): 
                    #print len(tvec), len(data_tvec) ,'len(tvec), len(data_tvec) '
                    all_data[:,los_num[i]] = interp(tvec,data_tvec[valid],data[:,i].data[valid])
                elif any(valid):
                    all_data[:,los_num[i]] = data[:,i].data
                    
                detector_stat[los_num[i]]  = self.status[sig]&~self.different_det[sig]

            detector_lim[los_num]   = ADCrange*calib+offset
            discr_err[los_num]      = calib
            detector_sampl[los_num] = SampFreq


            all_index.append(isig+arange(len(signals)))
            
        #import IPython
        #IPython.embed()
        #store data  
        #savez_compressed('tvec_'+str(self.shot), tvec=tvec)
        #for cam, ind in self.cam_ind.iteritems():
            #save(cam+"_"+str(self.shot), all_data[:,ind])
            #save(cam+'_%d_DASlim.npy'%self.shot, detector_lim[ind])
            #save(cam+'_%d_stat.npy'%self.shot, detector_stat[ind])
            #save(cam+'_%d_SampFreq.npy'%self.shot, detector_sampl[ind])
            #save(cam+'_%d_err.npy'%self.shot, discr_err[ind])

        #exit()
        
        
        self.dets_index = [array(v) for k,v in self.cam_ind.iteritems()]
        calib = ones(len(self.cam_ind))
        

        wrong_det = self.all_los[~detector_stat]
        self.wrong_dets_damaged = wrong_det

        #estimate errorbars!!
        nt = len(tvec)
        ind =  ~in1d(self.all_los,self.wrong_dets_damaged)
        wrong_data_ind = isnan(all_data)
        all_data[wrong_data_ind] = 0
        
  
        
        ind_t = slice(0,nt)
        if nt > 2000:  ind_t = r_[0, unique(random.randint(nt,size=1000)),nt-1]
        
        U,S,V = fast_svd(all_data[ind_t][:,ind],min(30,nt/3))

        
        
        #assume that the differnce between data and SVD retrofit is only noise 
        svd_err = ones(len(self.dets))*infty
        svd_err[ind] = std(all_data[ind_t][:,ind]-dot(U, V*S[:,None]),0)

        wrong_data_ind |= all(all_data==mean(all_data,1)[:,None],1)[:,None]
        wrong_data_ind |= all_data > detector_lim-1 #overburned

        #constant noise for all timepoints!
        all_data_err = zeros_like(all_data)
        all_data_err[:] = 2*svd_err+discr_err

        all_data_err[wrong_data_ind] = infty
    


        #mostly overburned values and missing detectors - just for plotting!
        for i in range(self.nl):  
            if any(isnan(all_data[:,i])):
                all_data[:,i][isnan(all_data[:,i])] = detector_lim[i]
   
        #apply special corrections
        self.hardcoded_corrections(tvec, all_data,all_data_err)


        return tvec, all_data, all_data_err


      
    
    
        

    def load_geom(self, path):
        
        
        #separate equlibrium for each campaign 
        ##load corrections in degrees of the camera position

        self.geometry_version = 0
        corrections = 'det_pos_corr_null'
        
        
        if self.shot> 34200:#jus guess
            self.geometry_version = 10
            corrections =  'det_pos_corr_2017_3'


        elif self.shot> 34007:
            self.geometry_version = 9
            corrections =  'det_pos_corr_2017_2'


        elif self.shot> 33800:
            self.geometry_version = 8
            corrections =  'det_pos_corr_2017'

        elif self.shot> 31776:
            self.geometry_version = 7
            corrections =  'det_pos_corr_new'

        elif self.shot> 27439:
            self.geometry_version = 6
            corrections =  'det_pos_corr'

        elif self.shot> 25600:
            self.geometry_version = 5
            corrections =  'det_pos_corr_old'

        elif self.shot> 24916:
            self.geometry_version = 4
            corrections =  'det_pos_corr_old3'

        elif self.shot>18000:
            self.geometry_version = 3
            
        elif self.shot> 13000:
            self.geometry_version = 2
            
        elif self.shot< 1000:  #artificial discharges!!
            self.geometry_version = 1
            corrections = 'det_pos_corr_null'

            
            


        pos_corr = loadtxt(path+'/'+corrections+'.txt',
                           dtype={'names': ('det', 'angle'),'formats': ('S2',  'f4')})
        pos_corr =  {k:item for k,item in pos_corr}

        Phi = []
        coord_dict = OrderedDict()
        
        #print pos_corr

        #cams = array(['F','G', 'H', 'I', 'J', 'K','L', 'M'])
        #f,axis = subplots(2,4,sharex=True,sharey=True)
        #f.subplots_adjust(hspace=0.10, wspace = 0.1)
        #from matplotlib.collections import PolyCollection
        #prepare files with cordinates
        for icam,(det,signals) in enumerate(self.detectors_dict.iteritems()):
            xfile = open(self.geometry_path+'/detector_%s_x.txt'%det,'w')
            yfile = open(self.geometry_path+'/detector_%s_y.txt'%det,'w')
            dist_file = open(self.geometry_path+'/dist_%s.txt'%det,'w')

            coord_dict[det] = []
            verts = []

            color = ('r', 'b', 'g', 'k', 'm', 'y', 'k', 'gray','r','m','k','b','g','y','k')
         
            m_alpha = []
            for sig in signals:
                r1 = self.geometry['RPINHOLE'][sig]
                z1 = self.geometry['ZPINHOLE'][sig]
                r2 = self.geometry['REND'][sig]
                z2 = self.geometry['ZEND'][sig]
                m_alpha.append(abs(abs(arctan2(z2-z1,r2-r1))-pi/2))
                Phi.append( self.geometry['Tor_Pos'][sig]+45)#add 45deg to make tor postion consistent with diaggeom
                
       
            m_alpha = mean(m_alpha)
         
         
            for sig in signals:
                r1 = self.geometry['RPINHOLE'][sig]
                z1 = self.geometry['ZPINHOLE'][sig]
                r2 = self.geometry['REND'][sig]
                z2 = self.geometry['ZEND'][sig]
                coord_dict[det].append([[r1,r2, z1,z2]])

                
                THETA = self.geometry['THETA'][sig]
                CAMANGLE  = self.geometry['CAMANGLE'][sig]
                
                dAngle = deg2rad(THETA-CAMANGLE)
                if dAngle > pi: dAngle-= 2*pi
                
     

                alpha = arctan2(z2-z1,r2-r1)
            
                if self.shot > 1000 and self.shot > 20000 and det  in pos_corr:
                    alpha+= deg2rad(pos_corr[det])
    

                L = hypot(r2-r1, z2-z1)
                if L == 0: L = 1;alpha = 1
                
     
                #camera geometry 
                Dx = self.geometry['D_Width'][sig] #  0.96 mm
                Delta = self.geometry['Foc_Len'][sig]#14.#mm
                Px = self.geometry['P_Width'][sig] # 0.3#mm
                Py = self.geometry['D_Length'][sig] # 4.6#mm
                Dy = self.geometry['P_Length'][sig] # 5.0#mm
                Delta1 = Px/Dx*Delta/(1+Px/Dx)
                theta = arctan(Px/2/Delta1)

                #angle between axis of the camera and LOS
                theta*= cos(dAngle)
                                
                delta = tan(Dy/2./Delta)
                delta*= 0.5
                dist = r1*sin(arctan(delta/cos(arctan((z2-z1+1e-6)/(r2-r1+1e-6)))))
                dist_file.write('%5.4f\n'%abs(dist))
                
                

                if m_alpha<pi/4:       
                    Lr = L*abs(sin(alpha)/sin(pi-theta-alpha))
                    Ll = L*abs(sin(pi-alpha)/sin(alpha-theta))

                    
                    xfile.write('%5.4f %5.4f %5.4f\n'%(r1+Ll*cos(alpha-theta),r1+Lr*cos(alpha+theta),r1))
                    yfile.write('%5.4f %5.4f\n'%(z1+Lr*sin(alpha+theta),z1))
                    
                    verts.append([[r1,z1],[r1+Ll*cos(alpha-theta),z1+Lr*sin(alpha+theta) ]
                                    ,[r1+Lr*cos(alpha+theta),z1+Lr*sin(alpha+theta)]])

                else:
                    if tan(pi-abs(alpha)-theta)*tan(pi-abs(alpha)+theta)<0 and sin(pi/2+alpha)>0:
                        #solve special case for almost exactly vertical LOS
                        z21 = z2  
                        z22 = z2+1e-2
                    else:
                        z21 = z1-(r2-r1)*tan(pi-abs(alpha)-theta)*sign(alpha)
                        z22 = z1-(r2-r1)*tan(pi-abs(alpha)+theta)*sign(alpha)
                     
                    yfile.write('%5.4f %5.4f %5.4f\n'%(z21,z22,z1))
                    xfile.write('%5.4f %5.4f\n'%(r2,r1))
                    verts.append([[r1,z1],[r2,z21],[r2,z22]])

            xfile.close()
            yfile.close()

            #verts = array(verts)
            #coll = PolyCollection(verts, color='b',edgecolor='b',alpha = .1)
            #ax = axis.flatten()[det[0] == cams][0]
            #ax.add_collection(coll)
            #ax.autoscale_view()
            #ax.axis('equal')
            #ax.set_xlim(1,2.2)
            #ax.set_ylim(-1.3,1.3)    
            #ax.text(1.93,-0.95,det[0])

        #show()

    
 

 

def main():
    import os
    import os,sys
    sys.path.append('/afs/ipp/home/t/todstrci/TRANSP/')
    
    sxr = loader_SXR(30382)
    T = time.time()
    sxr.load_geom()
    print time.time()-T

    tvec, data, error = sxr.get_signal_fast(4.6, 4.7)
    
    sxr.hardcoded_corrections(tvec, data)
    print 'finished'
 
    
    
    

if __name__ == "__main__":
    main()
