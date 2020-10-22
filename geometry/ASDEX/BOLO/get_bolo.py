#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
    

import os,sys
from collections import OrderedDict
from numpy import *
from matplotlib.pyplot import *
from annulus import get_bd_mat



def ReadTomo(path, diag, shot, time,geometry_path):

    f = open(path+'resultef2d.%s.%d_%.6f'%(diag, shot, time))
    diaf, shot_, time_ = f.readline().split()
    ed = f.readline()
    nx,ny = f.readline().split()
    nx,ny = int(nx)+1, int(ny)+1
    R = array([float(f.readline().replace('D','E')) for i in range(nx)])
    z = array([float(f.readline().replace('D','E')) for i in range(ny)])
    E = [float(f.readline()) for i in range(ny*nx)]
    dR = mean(diff(R))
    dz = mean(diff(z))

    E = reshape(E, (nx,ny)).T
    ndets  = int(f.readline())
    data,retro,err,fact = array([float_(f.readline().split()) for i in range(ndets)]).T
    err[err<0] = infty
    err*= data
    r1, r2 = f.readline().split()
    chi2 = float(f.readline())


    f = open(path+'paramsef2d.BLB.%d_%.6f'%(shot, time))

    nx,ny,ver = f.readline().split()
    nx,ny = int(nx)+1, int(ny)+1

    R_ = [float(f.readline()) for i in range(nx)]
    z_ = [float(f.readline()) for i in range(ny)]

    psi,par,per = array([float_(f.readline().split()) for i in range(ny*nx)]).T
    psi = reshape(psi, (nx,ny))
    par = reshape(par, (nx,ny))
    per = reshape(per, (nx,ny))

    ed2 = f.readline()
    ndets = int(f.readline())

    r1,z1,psi,d,theta,f,_,_,_,_,data_,gain,used,fact_ = array([float_(f.readline().split()) for i in range(ndets)]).T


    border = loadtxt(geometry_path+'/border.txt')
    from scipy.ndimage.interpolation import zoom

    #zoom to reduce discretisation error when the radiation outside of the borders is set to zero
    n = 10
    E = zoom(E, n,order=0)
    R = zoom(R, n,order=0)
    z = zoom(z, n,order=0)
    dR = mean(diff(R))
    dz = mean(diff(z))
    

    bbox = R[0]-dR/2, R[-1]+dR/2, z[0]-dz/2, z[-1]+dz/2
    BdMat = get_bd_mat(tokamak=None, nx=nx*n, ny=ny*n, boundary=border, bbox=bbox)
    
    E[BdMat.reshape(ny*n, nx*n,order='F')] = 0


    power = 2*pi*sum(E*R[None])*dR*dz



    return R, z, E, data_*fact,retro,err, power, bbox, border,fact




    
class loader_BOLO():

    
    """
    class for the loading of the foil bolometer data from AUG

    """
    
    
    dt = 0.001
    sigma = 0.00
    min_error = 0.02
    
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

        #dd = dd.shotfile()
       
        calib_shot = self.dd.cShotNr('AUGD', 'BLC',self.shot)
       
        if not self.dd.Open( 'BLC', calib_shot):
            #second try
            calib_shot = self.dd.cShotNr('AUGD', 'BLC', calib_shot-1)
            if not self.dd.Open( 'BLC', calib_shot):
                raise Exception('No BLC shotfile!')
            
       
       
       
    
       
        names = [n.strip() for n in self.dd.GetNames()]

        cams = [b for b in names if b[0] == 'F' and len(b) == 3]
   
        active_los = {s:bool_(self.dd.GetParameter(s, 'active' )) for s in cams}
        activ_cam = [c for c in cams if any(active_los[c])]
        try:
            activ_cam.remove('FEC')
        except:
            pass

        
        self.Phi     = {s:self.dd.GetParameter(s, 'P_Blende' )[active_los[s]] for s in activ_cam}
        self.R_start = {s:self.dd.GetParameter(s, 'R_Blende' )[active_los[s]] for s in activ_cam}
        self.z_start = {s:self.dd.GetParameter(s, 'z_Blende' )[active_los[s]] for s in activ_cam}
        self.R_end   = {s:self.dd.GetParameter(s, 'R_end'    )[active_los[s]] for s in activ_cam}
        self.z_end   = {s:self.dd.GetParameter(s, 'z_end'    )[active_los[s]] for s in activ_cam}
        self.delta   = {s:self.dd.GetParameter(s, 'delta'    )[active_los[s]] for s in activ_cam}
        self.nchan   = {s:self.dd.GetParameter(s, 'N_chan') for s in activ_cam}
        self.raw_name= {s:self.dd.GetParameter(s, 'RAW') for s in activ_cam}

        self.dd.Close()

        self.theta = {s:arctan2(self.z_end[s]-self.z_start[s],self.R_end[s]-self.R_start[s])  for s in activ_cam}
        self.subcam_ind = {c:[self.R_start[c] == r for r in unique(self.R_start[c])] for c in activ_cam}
        self.channels = {s:arange(len(a))[a]  for s,a in  active_los.items()}

     
        self.activ_cam = activ_cam
        self.active_los = active_los
   
        self.detectors_dict = OrderedDict()
        
        self.cam_ind = OrderedDict()
        nlos = 0
        self.all_los = []
        for c in activ_cam:
            self.cam_ind[c] = nlos+arange(sum(active_los[c]))
            nlos+= sum(active_los[c])
            self.detectors_dict[c] = [c+'_%d'%ch for ch in self.channels[c]]
            self.all_los.extend(self.detectors_dict[c])
        self.nl = nlos


        self.calb_0 = ones(len(activ_cam))
        self.dets = arange(self.nl)
        self.dets_index = [v for k,v in self.cam_ind.items()]
        
        
        if not self.dd.Open('BLB',self.shot, experiment=self.experiment, edition=self.ed):
            raise Exception('Bolometric data do not exists!!!')
        tvec = self.dd.GetTimebase('timFVC')
        self.tvec = arange(tvec[0],tvec[-1], self.dt )  #timevector for dowsampled data
        
        if self.shot > 33800 :
            self.wrong_dets_damaged += [  'FDC_4',]
        
    def get_total_rad(self):
        
        
        def reduce(x,y,n=1000):
            r = len(x)//n
            if r > 1:
                nt = len(x)
                x = mean(x[:(nt//r)*r].reshape(-1, r),1)
                y = mean(y[:(nt//r)*r].reshape(-1, r),1)

            return x,y

                
        out = []
    
        if self.dd.Open('TOT',self.shot):
            P_TOT = self.dd.GetSignal('P_TOT')
            P_TOT_tvec = self.dd.GetTimebase('P_TOT')
            self.dd.Close()
            P_TOT_tvec,P_TOT = reduce(P_TOT_tvec,P_TOT)
            out.append(['TOT:P_TOT', P_TOT_tvec,P_TOT])

        
        if  self.dd.Open('BPD',self.shot):
        
            for sig in ('Pradtot','Prad',):
                if sig in self.dd.GetNames():
                    Prad = self.dd.GetSignal(sig)
                    Prad_tvec = self.dd.GetTimebase(sig)
                    out.append(['BPD:'+sig, Prad_tvec,Prad])
    
            self.dd.Close()
            

        
        from os import listdir
        path = '/afs/ipp-garching.mpg.de/home/m/mbern/entfaltungen/data/'

        folders = listdir(path)
        folders.sort()
        shot_folders = {}
        ed  = 0
        times= []
        diag = None
        for folder in folders:
            try:
                diag_, shot_ = folder.split('.')
                shot_, ed_ = shot_.split('_')
                shot_, ed_ = int(shot_),int(ed_)
            except:
                continue
            

            if self.shot != shot_ :continue


            files =  listdir(path+'/'+folder)
            times_ = []
            for file in files:
                if file.startswith('resultef2d'):
                    times_.append(float(file.split('_')[1]))
            
            if len(times_) > len(times) or ((len(times_) == len(times)) and (ed_>ed)):
                times = times_
                ed = ed_
                diag = diag_
        
        if not diag is None:
            total_power = []
            for time in times:
                folder = path+diag+'.%d_%d/'%(self.shot,ed)
                R, z, E, data,retro,err, power,bbox, border,fact = ReadTomo(folder, diag, self.shot, time, self.geometry_path)
                total_power.append(power)
            print( "loaded M. Bernett's tomography ed:%d  times:"%ed, times)
            out.append(['BLB:mbern', times,array(total_power)])

                        
            #import IPython
            #IPython.embed()
        
            #used_los = where(hstack([ self.active_los[c][:self.nchan[c]] for c in self.activ_cam]))[0]
            #self.corrections = fact[used_los]

            #tvec, data_all, data_err_all = self.get_data(tmin=time-.1,tmax=time+.1)

            #f, ax = subplots(1,2, figsize=(12,6))
            #ax[0].plot(retro/1e6)

            #ax[0].errorbar(arange(len(data)),data/1e6,err/1e6 ,c='r' )
            #ax[0].set_ylabel('Brightness [MW/m$^2$]')
            #ax[0].set_ylim(0,nanmax(data)*1.2/1e6)
            
            
            #ax[0].plot(used_los,nanmean(data_all,0)/1e6,'k', lw=.3)  #BUG different correction factors!!! 

            #ax[1].set_title('power %.2f MW'%(power/1e6))
            #im = ax[1].imshow(E/1e6,origin='lower',vmin=0,interpolation='nearest', extent=bbox)
            #ax[1].contour(   R,z, E, levels=(-1e-3,),colors='w')
            #ax[1].plot(border[:,0], border[:,1],'w')
            #ax[1].axis(bbox)
            #colorbar(im)
            #tight_layout()
            #savefig(os.path.expanduser('~/tomography/tmp/mbern_rad_power.png'),bbox_inches='tight')
    
            #close()
      
                
        return out
        
        
      
      
    def get_data(self,tmin=-infty,tmax=infty):

        if not self.dd.Open('BLB',self.shot, experiment=self.experiment, edition=self.ed):
            raise Exception('no bolometric data!!')
        
        from scipy.interpolate import interp1d
        
        imin,imax = self.tvec.searchsorted((tmin, tmax))
        imax += 1
        tvec = self.tvec[imin:imax]
        ntim = len(tvec)

        
        data_all = empty((ntim, self.nl),dtype=single)
        data_err_all = empty((ntim, self.nl),dtype=single)

        for i,det in enumerate(self.activ_cam):
            data_tvec = self.dd.GetTimebase('tim'+det)
            nbeg,nend = data_tvec.searchsorted([tmin,tmax])
            data_tvec = data_tvec[nbeg:nend]
            data = self.dd.GetSignal('pow'+det,True, nbeg=nbeg+1,nend=nend)
            #print( det, data.shape, shape(self.active_los[det]), self.nchan[det])
            #print(self.active_los[det])
            data = data[:,:self.nchan[det]][:,self.active_los[det][:self.nchan[det]]]
            
            
            #raw_data = [self.dd.GetSignal('raw'+det,False, nbeg=nbeg+1,nend=nend)]

            
            #raw_data = self.dd.GetSignal('raw'+det,False, nbeg=nbeg+1,nend=nend)
            #raw_data = raw_data[:,self.active_los[det]]
            
            #ene = self.dd.GetSignal('ene'+det,True, nbeg=nbeg+1,nend=nend)
            #ene = ene[:,self.active_los[det]]
            
             #self.raw_name[det][31]
             
             #self.dd.Open('BLV', self.shot)
                #ene = self.dd.GetSignal(self.raw_name[det][33])


            
            
            
            #ene[~isfinite(ene)] = 2e6
            #for i in range(ene.shape[1]):
                #ene[:,i] = convolve(ene[:,i]>1e6, ones(min(100,nend-nbeg)),'same')
            #data[bool_(ene)] = nan  #corrupted points
            
            #plot( data[:,-1]);show()
            
            #import IPython
            #IPython.embed()
            
            #self.raw_name
            
            
            #print det, sum(isnan(data))

            #dowsample to 0.001s resolution
            nt = nend-nbeg
            reduct = int(ceil(self.dt/(tmax-tmin)*nt))
            nt = (nt//reduct)*reduct
            data = reshape(data[:nt,:],(nt//reduct,reduct,-1))
            dataerr = std(data,axis=1)
            data = mean(data,axis=1)
            
            data_tvec = reshape(data_tvec[:nt],(nt//reduct,reduct))
            data_tvec = data_tvec.mean(1)

            if len(data_tvec) > 1:
                data_tvec[[0,-1]] = tmin-1,tmax+1
                data = interp1d(data_tvec,data,fill_value =0,bounds_error=True,
                            copy=False,axis=0,assume_sorted=True)(tvec)
                dataerr = interp1d(data_tvec,dataerr,fill_value=0,bounds_error=True,
                            copy=False,axis=0,assume_sorted=True)(tvec)

            dataerr[isnan(data)] = infty
            data[isnan(data)] = 0
            data_all[:,  self.cam_ind[det]] = data
            data_err_all[:, self.cam_ind[det]] = dataerr
         

        self.dd.Close()
   
        def apply_correction(cam,channels,corr):
            if cam in self.cam_ind:
                ind =  in1d(self.channels[cam], channels)
                data_all[:,self.cam_ind[cam][ind]]*= corr
                data_err_all[:,self.cam_ind[cam][ind]]*= corr

        data_err_all+= nanmean(data_all,1)[:,None]*self.min_error+self.sigma*data_all
        #BUG adhoc corrections of the foil bolometers

        if self.shot > 31591: #BUG find the accurate upper boundary!!
            apply_correction( 'FHC',  r_[12:17], 1.3)
            if self.shot < 32000:
                apply_correction( 'FVC', r_[16:49], 0.58)
        if self.shot > 30135:
                
            apply_correction( 'FVC',  16, 2.4)
            apply_correction( 'FVC',  r_[0:48], 0.87)
            apply_correction( 'FHC',  r_[36:48], 0.76)
            apply_correction( 'FHC',  33, 0.91)
            apply_correction( 'FHC',  34, 0.95)
            apply_correction( 'FHC',  r_[20:24], 1.2)
            ##apply_correction( 'FHC',  r_[0:12],0.7)
            apply_correction( 'FHC',  r_[0:48], 0.87)
            apply_correction( 'FHC',  32, 1.1)
            apply_correction( 'FHC',  [22,25,30], 1.05)

                
        elif self.shot > 28523:
            apply_correction( 'FVC',  r_[16:48], 1/1.75)
            apply_correction( 'FHC',  r_[36:48], 1/1.1)
            apply_correction( 'FHC',  [33,36], 1/1.1)


        elif self.shot > 27352:
            apply_correction( 'FHC',  33, 1/1.2)
            apply_correction( 'FHC',  34, 1/1.1)

        apply_correction( 'FHS', r_[0:48],1.1)

        #if hasattr(self, 'corrections')
            #data_all*= self.corrections

        return  tvec, data_all, data_err_all

    
   
   
   
   
       
 
    def load_geom(self,path):
        
        if self.shot < 26000:  #just guess!!
            self.geometry_version = 4
        elif self.shot < 29100:  #just guess!!
            self.geometry_version = 3
        elif self.shot < 30162 : 
            self.geometry_version = 4
        elif self.shot < 31776 : 
            self.geometry_version = 2
        else: 
            self.geometry_version = 1

        for i,det in enumerate(self.activ_cam):

            R1 = self.R_start[det]
            Z1 = self.z_start[det]
            R2 = self.R_end[det] 
            Z2 = self.z_end[det]

            #move the begining of LOS to the center of the fan
            I = ones_like(R1)
            if det == 'FLX': R1,R2, Z1,Z2 = 1.899*I,R1,-0.945*I,Z1
            if det == 'FDO': R1,R2, Z1,Z2 = 1.473*I,R2,-1.173*I,Z2
            if det == 'FDI': R1,R2, Z1,Z2 = 1.313*I,minimum(R1,R2),-1.148*I,maximum(Z1,Z2)

            
            m_alpha = abs(abs(arctan2(Z2-Z1,R2-R1))-pi/2)
            alphas = arctan2(Z2-Z1,R2-R1)
            alphas = unwrap(alphas)
            m_alpha = mean(m_alpha)
            
            
            xfile = open(self.geometry_path+'/detector_%s_x.txt'%det,'w')
            yfile = open(self.geometry_path+'/detector_%s_y.txt'%det,'w')
            
            
            alpha = arctan2(Z2-Z1,R2-R1) #to same jako alphas? 
            
            theta = ones_like(alpha)*median(abs(diff(alphas)))/2 #exactly nonoverlaping LOS 
            theta/= 4  #actual LOS are narrower...
                        
            if det == 'FHC':
                theta[:9] = median(abs(diff(alphas))[:9])/2
                theta[9:32] = median(abs(diff(alphas))[9:32])/2
                theta[32:] = median(abs(diff(alphas))[32:])/2


            L = hypot(R2-R1, Z2-Z1)

            verts = []

            if m_alpha < pi/4:       
                Lr = L*abs(sin(alpha)/sin(pi-theta-alpha))
                Ll = L*abs(sin(pi-alpha)/sin(alpha-theta))
                R21 = R1+Ll*cos(alpha-theta)
                R22 = R1+Lr*cos(alpha+theta)
                
                savetxt(self.geometry_path+'/detector_%s_x.txt'%det, c_[R21,R22,R1], fmt='%5.4f')
                savetxt(self.geometry_path+'/detector_%s_y.txt'%det, c_[Z2,Z1], fmt='%5.4f')
                for r1,r21,z1,r22,z2 in zip(R1,R21,Z1,R22,Z2):
                    verts.append([[r1,z1],[r21,z2 ],[r22,z2]])

            else:
                
                
                #ind = tan(pi-abs(alpha)-theta)*tan(pi-abs(alpha)+theta) >= 0
                 ###solve special case for almost exactly vertical LOS
                #Z21 = copy(Z2)
                #Z22 = copy(Z2) +1e-2
                
                #Z21[ind] = (Z1-(R2-R1)*tan(pi-abs(alpha)-theta)*sign(alpha))[ind]
                #Z22[ind] = (Z1-(R2-R1)*tan(pi-abs(alpha)+theta)*sign(alpha))[ind]
                
                Z21 = (Z1-(R2-R1)*tan(pi-abs(alpha)-theta)*sign(alpha))
                Z22 = (Z1-(R2-R1)*tan(pi-abs(alpha)+theta)*sign(alpha))
                
                
                
                

                savetxt(self.geometry_path+'/detector_%s_x.txt'%det, c_[R2,R1], fmt='%5.4f')
                savetxt(self.geometry_path+'/detector_%s_y.txt'%det, c_[Z21,Z22,Z1], fmt='%5.4f')  
   
                for z1,z21,r1,z22,r2 in zip(Z1,Z21,R1,Z22,R2):
                    verts.append([[r1,z1],[r2,z21],[r2,z22]])


            #color = ('r', 'b', 'g', 'k', 'm', 'y', 'k', 'gray')

            #from matplotlib.collections import PolyCollection
            #title(det)
            #verts = array(verts)
            #ax = gca()
            #coll = PolyCollection(verts, color=color[i],edgecolors='none',alpha = 0.1)
            #ax.add_collection(coll)
            #ax.autoscale_view()
            #ax.axis('equal')
            #ax.set_xlim(1.1, 2.15)
            #show()
    


 

def main():
    import os
    import os,sys
    sys.path.append('/afs/ipp/home/t/todstrci/TRANSP/')
    
    bolo = loader_BOLO(30287)
    T = time.time()
    bolo.get_data(4,6)
    bolo.get_total_rad(4,6)

    print( 'loaded in ',time.time()-T)
    bolo.load_geom()

    #tvec, data, error = sxr.get_signal_fast(4.6, 4.7)
    
    #sxr.hardcoded_corrections(tvec, data)
    #print 'finished'
 
    #ind = argmax(data.max(0))
    #errorbar(tvec, data[:,ind], error[:,ind])
    #show()
    
    
    
    
    

if __name__ == "__main__":
    main()
