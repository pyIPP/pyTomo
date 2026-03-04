#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
import sys,os
from matplotlib.pyplot import *
from tokamak import Tokamak
from collections import OrderedDict
import config
#AUTOR: Tomas Odstrcil  tomas.odstrcil@ipp.mpg.de

#BUG xtomo ignuruje ten čas začátku a konce!!
 
class TCV(Tokamak):
    """
    Representation of class of tokamak. Contain basic information about geometry, 
    borde, magnetic field and basic setting used for recontruction and postprocessing.

    >>> TCV(input_path, 'XTOMO',  shot, nx, ny,virt_chord )


    :var str name:  Name of tokamak
    :var int index: Index of the tokamak
    :var double xmin, xmax, ymin, ymax: boundarys of recontructed area of tokamak emissivity
    :var str geometry_path: Path to saved geometry, boundary or data
    :var double sigma:  Expected noise in data  += ``sigma``* sqrt(data) * sqrt(max(data))  (Poisson like noise)
    :var double min_error: Expected noise in data  += ``min_error`` * sqrt(max(data))
    :var int shot:  NUmber of selected shot
    :var int norm: Used to normalize input parameters to meters
    :var str t_name:  Physical units used to measure time for the tokamak
    :var int npix:  total number of pixels (nx*ny)
    :var double dx: Width of one cell
    :var double dy: Height of one cell
    :var bool local_network: True if it is possible to load data from network

    :var array data: Input signals from main diagnostics
    :var array tvec: Time vector that belongs to the ``data``
    :var array tsteps: Lenght of tvec /  data
    :var array error: Expected errors of data,
    :var array dets: Array of all detectors in diagnostics
    :var array wrong_dets: Array of deleted detectors in diagnostics
    :var array magx, magy: Arrays of magnetic lines possitions
    :var spmatrix Tmat: Matrix of geometry
    :var array Xchords, Ychords: Arrays of chords coordinates, used in ``make_plots``
    :var array boundary: Coordinates of points on boundary of tokamak
    """


    def __init__(self, input_diagn, input_parameters,  only_prepare = False  ):

        """ Initialize new tokamak and set default values.

            :param dict input_parameters:   List of selected paramerers
            :param str input_diagn:  Name of main diagnostics
            :param array magfield_shift: artificial shift of magnetic field

        """
        import socket
        self.remote_run = True
        if socket.gethostname().find('lac') != -1:
            self.remote_run = False
        
        
        
        name  = 'TCV'

        self.mds_server = input_parameters['mds_server']
        self.rad_profile = 'Gaussian'

        print(input_diagn)
        tcv_diags = ['XTOMO','AXUV','BOLO','DMPX','FIR','XTEPRO']
        self.shot = input_parameters['shot']
        self.load_mag_equilibrium(input_diagn)

        if input_diagn in tcv_diags:
            coord= [0.624, 1.136, -0.75, 0.75]
            
            if input_diagn in ['DMPX','XTOMO','FIR','XTEPRO']:
                #smallest boundary contaning the whole plasma
                mag_boundary = [0.624, 1.136,round(self.magy.min()-.05,2),
                                round(self.magy.max()+.05,2)]
                coord = mag_boundary

        else:
            raise Exception('diag was not implemented yet')

        
        self.index = [i+15 for i,d in enumerate(tcv_diags) if  input_diagn==d][0]

        path = 'geometry/'+name+'/'+ input_diagn

        self.allow_negative = False
        sigma = 0.1
        min_error = 0.01
        norm = 1.0
        t_name = 's'
        
        
        self.min_tvec = 0.
        self.max_tvec = 2.6

       
        Tokamak.__init__(self, input_parameters, input_diagn, coord,
                         name, norm, sigma, min_error, path,t_name)


        #self.XTEPRO()
        
        #exit()

        if input_diagn == 'XTOMO':
            self.max_tvec = 2
            self.XTOMO()
        if input_diagn == 'BOLO':
            self.BOLO()     
        if input_diagn == 'DMPX':
            self.DMPX() 
        if input_diagn == 'AXUV':
            self.AXUV()
        if input_diagn == 'FIR':
            self.FIR()
        if input_diagn == 'XTEPRO':
            self.XTEPRO()
            
        self.VesselStructures()


        Tokamak._post_init(self)
        self.Tokamak = Tokamak


        self.get_injection_times()

        if not only_prepare:
            self.prepare_tokamak()
            
            

        print("object TCV done")
        
        self.input_diagn =input_diagn 
        
        
        

        
    def XTOMO(self):

        #BUG kalibrace - mám etendue, ale chybí mi konvezní faktor z V do A a pak ddo W a pak na plochu diody,, W/m^2


        
        n_diods = 20
        n_dets = 10
        self.detectors_dict=OrderedDict()
        
        if self.shot > 49000:
            self.calb_0 = [0.822, 0.902, 0.849, 0.898, 0.563, 0.858, 0.966, 1.022, 1.172, 1.216]
            
        else:
            self.calb_0 = ones(n_dets)
        self.nl =  n_diods*n_dets
        
        self.dets_index = [r_[n_diods*i:n_diods*(i+1)]for i in range(n_dets)]
        self.dets = arange(self.nl)
        
        #list of permantly wrong channels
        config.wrong_dets_pref = r_[config.wrong_dets_pref, 163 -1]
        
        
        #import IPython
        #IPython.embed()
        
        


        try:
            assert self.remote_run
            #pp
            data = load(self.geometry_path+'/data_%d.npz'%self.shot)
            tvec = data['tvec']
            saturated = data['saturated'] if 'saturated' in data else None 
            data = data['data']



        except :
            from .XTOMO import get_xtomo
            print('load xtomo')
            import MDSplus as mds
            c = mds.Connection(self.mds_server )
            #


            #pmds.mdsconnect( self.mds_server)
            print('connected to '+self.mds_server)

            tvec,data,saturated,geom = get_xtomo.dget_xtomo(c,self.shot,tmax=self.tsurf[-1]+.1)
            #pmds.mdsdisconnect()
            c.closeAllTrees()

            #import IPython
            #IPython.embed()

            if self.remote_run:
                savez(self.geometry_path+'/data_%d.npz'%self.shot,tvec=tvec,
                    data=data,dets=self.dets,saturated=saturated)
            geom_dict = {}

            for idet, ind in enumerate(self.dets_index):
                cam_name = chr(ord('A')+idet)       
                xchord = geom['xchord'][:,ind].T/100
                ychord = geom['ychord'][:,ind].T/100
                self.LOSextend(cam_name, xchord, ychord)
                geom_dict[cam_name] = (xchord, ychord)
                
            save(self.geometry_path+'/geom_%d.npy'%self.shot, geom_dict)


        
        for idet, ind in enumerate(self.dets_index):
            cam_name = chr(ord('A')+idet)        
            self.detectors_dict[cam_name] = [cam_name+'_%.2d'%i for i in range(n_diods)]



        self.tvec = tvec
    
        #BUG camera I is incorrecly wired!!
        data[:,160:180] = copy(data[:,179:159:-1])
        if not saturated is None:
            saturated[:,160:180] = copy(saturated[:,179:159:-1])
        ind = r_[arange(200)[::20],arange(200)[19::20]]
        data[:,ind]/= 1.1
        
        #data*=1e4  #probably convert to Wm^-2??

     

        
        data_err = data[:tvec.searchsorted(0)].std(0)
        print('saving')
        self.save_data(data.T,tvec=tvec,dets=self.dets,err=data_err,
                        saturated=transpose(saturated))



    def DMPX(self):

        self.detectors_dict=OrderedDict()
        
        sigma = 0.05
        min_error = 0.1
        
        n_dets = [64,32]
        n_cams = len(n_dets)
        self.nl = 0
        self.dets_index = []
        for nd in n_dets:
            self.dets_index.append(r_[self.nl:self.nl+nd])
            self.nl+= nd
            
        cam_name = ['TOP','BOT']
        self.calb_0 = ones(n_cams)
        self.dets = arange(self.nl)
        #config.wrong_dets_pref = r_[64:97]

        try:
            assert self.remote_run

            data = load(self.geometry_path+'/data_%d.npz'%self.shot)
            tvec = data['tvec']
            self.dets = data['dets']
            saturated = data['saturated'] if 'saturated' in data else None 
            data = data['data']
        except:
            #import MDSplus as mds

            from .DMPX import get_dmpx
            
            
            import MDSplus as mds
            c = mds.Connection(self.mds_server )
            print('connected to '+self.mds_server)

            #pmds.mdsconnect( self.mds_server)
            trange = [-.05, self.tsurf[-1]]
   

            dmpx = get_dmpx.mpxdata(c,self.shot,'sfgev',time=trange)
            #pmds.mdsdisconnect()
            #c.closeAllTrees()


            effectivity = {}
            effectivity['top'] = dmpx.top.eff_ev, dmpx.top.eff_response
            effectivity['bot'] = dmpx.bot.eff_ev, dmpx.bot.eff_response
            #detector
            
            clf()
            semilogx(dmpx.top.eff_ev, dmpx.top.eff_response,label='TOP '+dmpx.top.detector)
            semilogx(dmpx.bot.eff_ev, dmpx.bot.eff_response,label='BOTTOM '+dmpx.bot.detector)
            legend(loc='upper right')
            ylim(0,1)
            xlim(1e3,5e4)
            xlabel('E [eV]')
            ylabel('absorbtion efficiency')
            savefig(self.geometry_path+'/absorbtion_efficiency%d.pdf'%self.shot)
            clf()

            tvec = dmpx.tvec
            data = dmpx.signal
            saturated = dmpx.saturated
            
            #downsample 4x  to 20kHz
            nt = (len(tvec)/4)*4
            tvec = mean(tvec[:nt].reshape(nt/4,4),1)
            data = mean(data[:nt].reshape(nt/4,4,-1),1)
            saturated = any(saturated[:nt].reshape(nt/4,4,-1),1)


            
        
            zmax = 0.75
            dz = (dmpx.top.geom['slitZ']-dmpx.top.geom['wiresZ'])
            dr = (dmpx.top.geom['slitR']-dmpx.top.geom['wiresR'])

            t = (zmax-dmpx.top.geom['wiresZ'])/dz
            z_end_top = dmpx.top.geom['wiresZ']+t*dz
            R_end_top = dmpx.top.geom['wiresR']+t*dr

            z_start_top = dmpx.top.geom['wiresZ']
            R_start_top = dmpx.top.geom['wiresR']


            dz = (dmpx.bot.geom['slitZ']-dmpx.bot.geom['wiresZ'])
            dr = (dmpx.bot.geom['slitR']-dmpx.bot.geom['wiresR'])

            t = (zmax-dmpx.bot.geom['wiresZ'])/dz
            z_end_bot = dmpx.bot.geom['wiresZ']+t*dz
            R_end_bot = dmpx.bot.geom['wiresR']+t*dr

            z_start_bot = dmpx.bot.geom['wiresZ']
            R_start_bot = dmpx.bot.geom['wiresR']

            geom_dict = {}
            geom_dict['TOP'] = (c_[R_start_top,R_end_top],c_[z_start_top,z_end_top])
            geom_dict['BOT'] = (c_[R_start_bot,R_end_bot],c_[z_start_bot,z_end_bot])
            save(self.geometry_path+'/geom_%d.npy'%self.shot, geom_dict)
            if self.remote_run:

                savez(self.geometry_path+'/data_%d.npz'%self.shot,tvec=tvec,data=data,
                    saturated = saturated, dets=self.dets, effectivity=effectivity)

       
       
       
       
        self.detectors_dict['TOP'] = ['T_%.2d'%i for i in range(n_dets[0])]
        self.detectors_dict['BOT'] = ['B_%.2d'%i for i in range(n_dets[1])]


        self.load_geom()



        self.tvec = tvec
     
        data_err = data[:tvec.searchsorted(1e-2)].std(0)
        print('saving')
        self.save_data(data.T,tvec=tvec,dets=self.dets,err=data_err,
                       saturated=transpose(saturated))



        
        
                
    def FIR(self):
        #WARNING !!
        #not a proper density deconvolution program
        #monotonisity contrain not included!
        
        


        self.detectors_dict=OrderedDict()
        
        
        
        self.sigma = 0.05
        self.min_error = 0.01

        self.nl = 14
        self.dets_index = [arange(self.nl),]

        #mean relatice calibration for beginning of campaing 2015
        self.calb_0 = ones(self.nl)

      
        try:
            assert self.remote_run

            data = load(self.geometry_path+'/data_%d.npz'%self.shot)
            tvec = data['tvec']
            self.dets = data['dets']
            data = data['data']

        except:
            from .FIR import get_fir
            import MDSplus as mds
            
            #import MDSplus as mds
            c= mds.Connection(self.mds_server )
            print('connected to '+self.mds_server)

            tvec,data,R,Z = get_fir.get_fir(c,self.shot)

            #c.closeAllTrees()



        

            #pmds.mdsconnect( self.mds_server)


            #pmds.mdsdisconnect()


            self.dets = arange(self.nl)
            if self.remote_run:

                savez(self.geometry_path+'/data_%d.npz'%self.shot,tvec=tvec,data=data,dets=self.dets)

            savetxt(self.geometry_path+'/detector_FIR_x.txt',c_[R,R],fmt='%.4f')
            savetxt(self.geometry_path+'/detector_FIR_y.txt',tile(Z,(self.nl,1)),fmt='%.4f')
                

        

        
        self.detectors_dict['FIR'] = ['%.2d'%i for i in range( self.nl)]



        self.tvec = tvec
     
        print('saving')
        self.save_data(data,tvec=tvec,dets=self.dets)

        
        
          
        
    def AXUV(self):
        from .AXUV import get_axuv
        #import MDSplus as mds
        
        #lapší odhad errorbarů? relativní komponenty vysokou? 
        
        self.detectors_dict=OrderedDict()
        self.boundary_lcfs = False

        self.sigma = 0.05
        self.min_error = 0.01
        n_cams = 7
        n_diods = 20
        
        self.nl = n_cams*n_diods
        #n_dets
                
        self.calb_0 = ones(n_cams)
        #self.nl =  n_diods* n_cams
        
        self.dets_index = [r_[n_diods*i:n_diods*(i+1)]for i in range(n_cams)]
        self.dets = arange(self.nl)
        
        #print self.nl
        #exit()


        try:
            assert self.remote_run
            #pp
            print('loading')
            #import IPython
            #IPython.embed()
            data = load(self.geometry_path+'/data_%d.npz'%self.shot)
            tvec = data['tvec']
            self.dets = data['dets']
            saturated = data['saturated']
            data = data['data']
            self.load_geom()


        except:
            
            import MDSplus as mds
            c= mds.Connection(self.mds_server )
            print('connected to '+self.mds_server)
            bolostruct = get_axuv.axuv_dld(c,self.shot)

            #c.closeAllTrees()






            #pmds.mdsconnect( self.mds_server)
            
            tvec = bolostruct.time
            data = bolostruct.data
            saturated = bolostruct.saturated

            #pmds.mdsdisconnect()
            #nl = data.shape[1]
           
            #self.dets = array(nl)
            
        
            if self.remote_run:

                savez(self.geometry_path+'/data_%d.npz'%self.shot,tvec=tvec,data=data,
                      dets=self.dets,saturated=saturated)

            #for idet, ind in enumerate(self.dets_index):
                #savetxt(self.geometry_path+'/detector_'+str(idet)+'_x.txt',
                        #bolostruct.xchord[:,ind].T,fmt='%.4f')
                #savetxt(self.geometry_path+'/detector_'+str(idet)+'_y.txt',
                        #bolostruct.ychord[:,ind].T,fmt='%.4f')

            geom_dict = {}
            for idet, ind in enumerate(self.dets_index):
                cam_name = chr(ord('A')+idet)       
                xchord = bolostruct.xchord[:,ind].T
                ychord = bolostruct.ychord[:,ind].T
                self.LOSextend(cam_name, xchord, ychord)
                geom_dict[cam_name] = (xchord, ychord)
                
                save(self.geometry_path+'/geom_%d.npy'%self.shot, geom_dict)


  
          
        #for i in arange(n_dets):
            #self.detectors_dict['cam %.2d'%i] = ['cam_%.2d'%i for i in range(num)]

        for idet, ind in enumerate(self.dets_index):
            cam_name = chr(ord('A')+idet)       
            self.detectors_dict[cam_name] = ['cam %d'%(idet)+'_%.2d'%i for i in range(n_diods)]


        self.tvec = tvec
     
        data_err = data[:tvec.searchsorted(0)].std(0)
        print('saving')
        self.save_data(data.T,tvec=tvec,dets=self.dets,err=data_err,saturated=transpose(saturated))

        
        
    def XTEPRO(self):
        from .XTEPRO import get_xtepro
        
        self.detectors_dict=OrderedDict()

        sigma = 0.05
        min_error = 0.1
        
        n_dets = [16,16,16]
        n_cams = len(n_dets)
        self.nl = sum(n_dets)
        self.dets_index = []
        for i,nd in enumerate(n_dets):
            self.dets_index.append(r_[i*nd:nd*(i+1)])
            
        cam_name = ['low','med','hig']
        self.calb_0 = ones(n_cams)


        try:
            #pp
            assert self.remote_run

            data = load(self.geometry_path+'/data_%d.npz'%self.shot)
            tvec = data['tvec']
            self.dets = data['dets']
            data = data['data']

        except  Exception as e:

            
            import MDSplus as mds
            c= mds.Connection(self.mds_server )
            print('connected to '+self.mds_server)

            tvec,data,saturated,geom = get_xtepro.dget_xtepro(c,self.shot,tmax=self.tsurf[-1]+0.1)

            
            #c.closeAllTrees()

            self.dets = arange(self.nl)
            
            if self.remote_run:

                savez(self.geometry_path+'/data_%d.npz'%self.shot,tvec=tvec,data=data,dets=self.dets)

            #for idet, ind in enumerate(self.dets_index):
                #savetxt(self.geometry_path+'/detector_'+cam_name[idet]+'_x.txt',
                        #geom['xchord'][:,::-1].T,fmt='%.4f')
                #savetxt(self.geometry_path+'/detector_'+cam_name[idet]+'_y.txt',
                        #geom['ychord'][:,::-1].T,fmt='%.4f')
                
                            
            #import IPython
            #IPython.embed()
            geom_dict = {}
            xchord = geom['xchord'][:,::-1].T
            ychord = geom['ychord'][:,::-1].T
            for idet, name in enumerate(cam_name):
                self.LOSextend(cam_name, xchord, ychord)
                geom_dict[name] = (xchord, ychord)
                
            save(self.geometry_path+'/geom.npy', geom_dict)

        

        
        for name,num in zip(cam_name,n_dets):
            self.detectors_dict[name] = [name+'_%.2d'%i for i in range(num)]



        self.tvec = tvec
     
        data_err = data[:tvec.searchsorted(0)].std(0)
        print('saving')
        self.save_data(data.T,tvec=tvec,dets=self.dets,err=data_err)

        
        
        
        
    def BOLO(self):
        from .BOLO import get_bolo
        #import MDSplus as mds
        
        self.detectors_dict=OrderedDict()
        self.boundary_lcfs = False


        sigma = 0.05
        min_error = 0.1
        
        n_dets = [8,16,16,16,8]
        n_cams = len(n_dets)
        self.nl = 0
        self.dets_index = []
        for nd in n_dets:
            self.dets_index.append(r_[self.nl:self.nl+nd])
            self.nl+= nd
            
        cam_name = ['TV','TL','ML','BL','BV']
        self.calb_0 = ones(n_cams)


        try:
            assert self.remote_run

            data = load(self.geometry_path+'/data_%d.npz'%self.shot)
            tvec = data['tvec']
            self.dets = data['dets']
            data = data['data']

        except  Exception as e:

            
            import MDSplus as mds
            c= mds.Connection(self.mds_server )
            print('connected to '+self.mds_server)

            bolostruct = get_bolo.get_bolo(c,self.shot)


            #c.closeAllTrees()

            



            #pmds.mdsconnect( self.mds_server)
            #bolostruct = get_bolo.get_bolo(self.shot)
            
            tvec = bolostruct.time
            data = bolostruct.data

            #pmds.mdsdisconnect()
        

            self.dets = arange(self.nl)
            
            if self.remote_run:

                savez(self.geometry_path+'/data_%d.npz'%self.shot,tvec=tvec,data=data,dets=self.dets)

            for idet, ind in enumerate(self.dets_index):
                savetxt(self.geometry_path+'/detector_'+cam_name[idet]+'_x.txt',
                        bolostruct.geometry.xchord[:,ind].T,fmt='%.4f')
                savetxt(self.geometry_path+'/detector_'+cam_name[idet]+'_y.txt',
                        bolostruct.geometry.ychord[:,ind].T,fmt='%.4f')
                

        

        
        for name,num in zip(cam_name,n_dets):
            self.detectors_dict[name] = [name+'_%.2d'%i for i in range(num)]



        self.tvec = tvec
     
        data_err = data[:tvec.searchsorted(0)].std(0)
        print('saving')
        self.save_data(data.T,tvec=tvec,dets=self.dets,err=data_err)

        
        
        
        
        

    def ShotOverview(self):
        pass
        
        
    def VesselStructures(self):
        self.struct_dict = {}
        
        
        #inner boundary
        x = [1.137, 0.9685,0.678,0.624,0.624 ,0.667 ,0.9704 ,1.134,1.137 ]
        y = [0.547, 0.748,0.749,0.703,-0.704,-0.745,-0.747,-0.551,0.547]
      
        self.struct_dict['inner'] = x,y


        x = [1.154, 0.997,.6385,.6035,0.6,0.654, 1.012,1.154,1.154  ]
        y = [0.631,0.77, 0.769, 0.7241,-0.709,-0.77,-0.766,-0.631,0.631]
        
        self.struct_dict['outer'] = x,y

  
    def mag_equilibrium(self, tvec, **kwarg):
        #BUG very ugly solution :( should not be stored at all 
        rho,magx,magy = self.Tokamak.mag_equilibrium(self,tvec,**kwarg)


        return rho,magx,magy
        
        
        

    def load_mag_equilibrium(self,input_diagn):
        #  LOAD FIELD FROM NETWORK

        print("load_mag_equilibrium, ",self.shot)

        try:
            #pp
            d = load('geometry/TCV/equilibrium/MagField_%d.npz'%self.shot)

            self.tsurf = d['tsurf']
            self.magx = d['magx']
            self.magy = d['magy']
            
        except Exception as e:

            import MDSplus as mds
            c= mds.Connection(self.mds_server )

            print('loading equilibrium')
            from .get_psi import get_psi
            liuqenumber = 1  #BUG which one os the best?? 
            tsurf,mag_axis,lcfs,rho,cnt = get_psi(c,self.shot,liuqenumber,True)

            c.closeAllTrees()
                

            #pmds.mdsconnect(self.mds_server)  #adress of the MDS+ server

            #pmds.mdsdisconnect()  
            n_theta = 100
            n_rho = len(cnt[0])+1
 
                      
            self.magx = empty((n_theta, n_rho,len(tsurf)),dtype='single')
            self.magy = empty((n_theta, n_rho,len(tsurf)),dtype='single')
            self.tsurf = tsurf
            
            theta = linspace(-pi,pi,n_theta)
            self.magx[:,0] = mag_axis[0]
            self.magy[:,0] = mag_axis[1]
            
            
                    
        
            #import IPython
            #IPython.embed()
            
        


            
            
            #interpolate data on the same value of theta values
            for it, ct in enumerate(cnt):
                #print it
                for ir, cr in enumerate(ct[::-1]):
                    if len(cr) > 3:
                        theta_surf = arctan2(cr[:,0]-mag_axis[0][it],cr[:,1]-mag_axis[1][it])
                        ind = argsort(theta_surf)
                        cr_,theta_surf = cr[ind],theta_surf[ind]
                        theta_surf = r_[theta_surf[-1]-2*pi,theta_surf,theta_surf[0]+2*pi] #periodicity
                        cr_ = vstack((cr_[-1,],cr_,cr_[0,]))
                        self.magx[:,-ir-1,it] = interp(theta,theta_surf,cr_[:,0])
                        self.magy[:,-ir-1,it] = interp(theta,theta_surf,cr_[:,1])
                        #import IPython
                        #IPython.embed()
                        #plot(cr[:,0], cr[:,1],'+')
                        #bd = self.get_boundary(N=n_theta, coords =cr ,k=2)
                        ##plot(bd[:,0], bd[:,1])
                        #self.magx[:,-ir-1,it] = bd[:,0]
                        #self.magy[:,-ir-1,it] = bd[:,1]

                        

                    else:
                        #extrapolate this surface from the larger one
                        self.magx[:,-ir-1,it] = (self.magx[:,-ir,it]-mag_axis[0][it])/(n_rho-ir)*(n_rho-ir-1.)+mag_axis[0][it]
                        self.magy[:,-ir-1,it] = (self.magy[:,-ir,it]-mag_axis[1][it])/(n_rho-ir)*(n_rho-ir-1.)+mag_axis[1][it]

                show()
            savez_compressed('geometry/TCV/equilibrium/MagField_%d.npz'%self.shot, 
                    tsurf=self.tsurf, magx = self.magx, magy = self.magy)
            
            
        if input_diagn in ['XTOMO','DMPX','FIR']: 
       
            
            from shared_modules import inside_convex_hull
            from scipy.spatial import ConvexHull
            points = c_[self.magx[:,-1,:].flatten(),self.magy[:,-1,:].flatten()]
            hull = ConvexHull(points)
           
            separatrix = c_[points[hull.vertices,0], points[hull.vertices,1]]
            separatrix = vstack([separatrix[-1],separatrix])
            
            
            #from shared_modules import read_config
            #input_parameters = read_config('tomography'+".cfg")
            #magfield_shift = input_parameters['magfield_shift_lcfs']
            #if not input_parameters['only_core_shift']:
            #separatrix += magfield_shift

       
            #savetxt('geometry/TCV/'+ input_diagn+'/border.txt',separatrix,fmt='%.4e')
            #self.vessel_boundary = separatrix




    def LOSextend(self,cam,xchords, ychords,show_los=False):
        #only the estimate!! not base on the camera geometry!!
        
        ang_corr = None
        
        path = None
        if os.path.isfile(self.geometry_path +'/geom_corrections_%d.txt'%self.shot):
            path = self.geometry_path +'/geom_corrections_%d.txt'%self.shot
        elif os.path.isfile(self.geometry_path +'/geom_corrections.txt'):
            path = self.geometry_path +'/geom_corrections.txt'
        #path = None
        if not path is None:
            ang_corr = loadtxt(path,dtype={'names': ('det', 'angle'),'formats': ('S3',  'f4')})
            ang_corr =  {k:item for k,item in ang_corr}
        
        #print ang_corr
        #exit()
        

        from matplotlib.collections import PolyCollection

        verts = []
        #color = ('r','b','g','y','k','m')
        m_alpha = []
        alphas = []
        for [r1,r2],[z1,z2] in zip(xchords, ychords):
            m_alpha.append(abs(abs(arctan2(z2-z1,r2-r1))-pi/2))
            alphas.append(arctan2(z2-z1,r2-r1))

        alphas = unwrap(alphas)
        m_alpha = mean(m_alpha)
        

        xfile = open(self.geometry_path+'/detector_%s_x.txt'%cam,'w')
        yfile = open(self.geometry_path+'/detector_%s_y.txt'%cam,'w')
        for j,([r1,r2],[z1,z2]) in enumerate(zip(xchords, ychords)):
            alpha = arctan2(z2-z1,r2-r1)
            if show_los:text(r2,z2,j)

            if not ang_corr is None:
                #print cam, ang_corr[cam]/180*pi 
                alpha+= ang_corr[cam]/180*pi 

            theta = median(abs(diff(alphas)))/2 #exactly nonoverlaping LOS 

            L = hypot(r2-r1, z2-z1)


            if m_alpha < pi/4:       
                Lr = L*abs(sin(alpha)/sin(pi-theta-alpha))
                Ll = L*abs(sin(pi-alpha)/sin(alpha-theta))
                r21 = r1+Ll*cos(alpha-theta)
                r22 = r1+Lr*cos(alpha+theta)
                
                
                xfile.write('%5.4f %5.4f %5.4f\n'%(r21,r22,r1))
                yfile.write('%5.4f %5.4f\n'%(z2,z1))
                
                verts.append([[r1,z1],[r21,z2 ],[r22,z2]])

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
        if show_los:
            verts = array(verts)
            ax = gca()
            coll = PolyCollection(verts, color='r',edgecolors='none',alpha = 0.1)
            ax.add_collection(coll)
            ax.autoscale_view()
            title(self.input_diagn+' '+cam)
            coord= [0.624, 1.136, -0.75, 0.75]
            plot([0.624,1.136,1.136,0.624,0.624],[-0.75,-0.75,0.75,0.75,-0.75],'k')
            #ax.axis([0.624, 1.136, -0.75, 0.75])
            ax.axis('equal')
            show()




    def load_geom(self):
        #refresh the geometry (include new position correction)
        if os.path.isfile(self.geometry_path+'/geom.npy'):
            geom_dict = load(self.geometry_path+'/geom.npy').item()

        elif os.path.isfile(self.geometry_path+'/geom_%d.npy'%self.shot):
            geom_dict = load(self.geometry_path+'/geom_%d.npy'%self.shot).item()
        else:
            return 
            
        for k,g in geom_dict.items():
            self.LOSextend(k,g[0], g[1])
            
                
    def get_injection_times(self):
        
        
        
        self.impur_inject_t = None
        #self.impur_inject_t = loadtxt('geometry/TCV/impur_injections/%d.txt'%self.shot)

        try:
            self.impur_inject_t = loadtxt('geometry/TCV/impur_injections/%d.txt'%self.shot)
            self.impur_inject_t = array(self.impur_inject_t,ndmin=1)

        
        except:
            print('get_injection_times')
            import MDSplus as mds
            c= mds.Connection(self.mds_server )
            
            c.openTree('tcv_shot',self.shot)

            try:
                self.impur_inject_t = [nan]
                for i in ['12','34','56','78','9A','BC']:
                    ttxt=r'\vsystem::timers_diagdb_r["TIMER_86%s:ABSOLUTE"]'
                    t_on=c.get(ttxt % i[0])
                    t_off=c.get(ttxt % i[1])

                    if t_on>=0 and t_off >=0:
                        #times.append((t_on,t_off))
                        self.impur_inject_t.append(t_on)
                c.closeAllTrees()
                savetxt('geometry/TCV/impur_injections/%d.txt'%self.shot, self.impur_inject_t,fmt='%.5g')
            except:
                pass
        #return times

            
                
    def get_sigma_mat(self,time,dets):
        #estimate a correlation matrix of the data
    
        try:
            C = load(self.geometry_path+'/corr_matrix.npy')
            
            C = C[dets][:,dets]
            
            
            sigma = linalg.cholesky(C+diag(diag(C)))  #regularizovanou matici sigma??? 

            return  sigma
        
        
        except:
            return 
    

        
