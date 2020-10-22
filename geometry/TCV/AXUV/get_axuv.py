#!/usr/bin/env python
# -*- coding: utf-8 -*-


from numpy import *
#import pmds
#from matplotlib.pylab import *
#import ././tqdm
import os.path 

import sys,os
#print os.path.abspath('../../')
#sys.path.append(os.path.abspath('../../'))
#import tqdm
from scipy.io import loadmat


parent_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
print(parent_path)
sys.path.append(parent_path)
import tqdm


"""
python version of axuv_get_axuv.m

rewriten by Tomas Odstrcil   tom@cbox.cz




"""



class axuvstruct_class:
    def __init__(self):
        class geometry:
            pass
        self.geometry = geometry()


def get_axuv(connect, shotno,tbegin=None,tend=None, shotno_offset= None):
        
    axuvstruct = axuvstruct_class()
        
    bolomap = empty((20,7),dtype=int)
    bolomap[:,0]=[ 1,20, 2,19, 3,18, 4,17, 5,16, 6,15, 7,14, 8,13, 9,12,10,11]
    bolomap[:,1]=[43,62,44,61,45,60,46,59,47,58,48,57,49,56,50,55,51,54,52,53]
    bolomap[:,2]=[ 1,94, 2,93, 3,92, 4,91 ,5,90, 6,89, 7,88, 8,87, 9,86,10,85]
    bolomap[:,3]=[43,42,44,41,45,40,46,39,47,38,48,37,49,36,50,35,51,34,52,33]
    bolomap[:,4]=[85,84,86,83,87,82,88,81,89,80,90,79,91,78,92,77,93,76,94,75]
    bolomap[:,5]=[21,42,22,41,23,40,24,39,25,38,26,37,27,36,28,35,29,34,30,33]
    bolomap[:,6]=[65,84,66,83,67,82,68,81,69,80,70,79,71,78,72,77,73,76,74,75]
    boards = empty((20, 7),dtype=int)
    boards[:,r_[0:2, 3:7]]= [11, 11, 12, 12 ,13 ,13]
    boards[:,2]=[12,11,12,11,12,11,12,11,12,11,12,11,12,11,12,11,12,11,12,11]

    if shotno>=44499:  #the cabling changed after the recommissioning!!!
        bolomap[:,2:5] = bolomap[::-1,2:5]
        boards[:,2:5] = boards[::-1,2:5]
        axuvstruct.geometry=axuv_calc_los(axuv_detectorpos('axuv'))
        if os.path.isfile(diag_path+'axuv_finiteetend_axuv.mat'):
            geom = loadmat(diag_path+'axuv_finiteetend_axuv.mat')
            axuvstruct.geometry.etend=squeeze(geom['etendstruct']['etend'].item())
        
        gain=400e3
    elif shotno<44499 and shotno>=34666:
        axuvstruct.geometry=axuv_calc_los(axuv_detectorpos('axuv_old'))
        gain=20e3
    else:
        raise Exception('The routine is not suitable for downloading these data.')

    connect.openTree('tcv_shot',shotno)
    t0=float(connect.get('\ATLAS::DT196_AXUV_001:STOP_TRIG'))        #triggering time
    per=float(connect.get('slope_of(\ATLAS::DT196_AXUV_001:EXT_CLOCK_IN)'))  #acquisition period (maybe sampling time)
    ns=int(connect.get('\ATLAS::DT196_AXUV_001:MAX_SAMPLES*1000') )        #ns is number of samples acquired
    connect.closeTree('tcv_shot',shotno)
    #pmds.mdsopen('tcv_shot',shotno)
    #t0=pmds.mdsvalue('\ATLAS::DT196_AXUV_001:STOP_TRIG')        #triggering time
    #per=pmds.mdsvalue('slope_of(\ATLAS::DT196_AXUV_001:EXT_CLOCK_IN)')  #acquisition period (maybe sampling time)
    #ns=pmds.mdsvalue('\ATLAS::DT196_AXUV_001:MAX_SAMPLES*1000')         #ns is number of samples acquired
    #pmds.mdsclose('tcv_shot',shotno)
    
    
    
    time_full= t0+arange(ns)*per
    if time_full[0]<0:
        tminoffsetindex = 0
        tmaxoffsetindex = time_full.searchsorted(0)
    elif time_full[-1]>=2:
        tminoffsetindex= time_full.searchsorted(2)
        tmaxoffsetindex= len(time_full)-1
    else:
        tmaxoffsetindex=None
        tminoffsetindex=None


    if tbegin is None:
        tbegin=time_full[0]
        tbeginindex=None

    elif time_full[0]>tbegin:
        print('Data acquisition started at ' +str(time_full[0])+ ' sec.')
        tbegin=time_full[0]
        tbeginindex=None
    else:
        tbeginindex = time_full.searchsorted(tbegin)

    if tend is None:
        tend=time_full[-1]
        tendindex = None
    elif time_full[-1]<tend:
        print('Data acquisition finished at '+ str(time_full[-1])+' sec.')
        t_end=time_full[-1]
        tendindex = None
    else:
        tendindex = time_full.searchsorted(tend)+1
        
    t_ind = slice(tbeginindex,tendindex)
    if (tbeginindex is None and not tendindex is None):
        tbeginindex = 0
    if (not tbeginindex is None and tendindex is None):
        tendindex = len(time_full)
    #tendindex=t_ind.stop
    time=time_full[t_ind]


    data_raw=zeros((len(time),140),dtype='single')
    offset = zeros(140)
    #pmds.mdsopen('crpppc250.epfl.ch::axuv,$',shotno)
    #pmds.mdsopen('crpppc250.epfl.ch',shotno)
    #pmds.mdsopen('tcv_shot',shotno)
    connect.openTree('tcv_shot',shotno)

    #if tmaxoffsetindex is None:
    for i in range(data_raw.shape[1]):
        data_raw[:,i],offset[i]=get_chord_offset(connect,boards.T.flat[i],bolomap.T.flat[i],
            tbeginindex,tendindex,tminoffsetindex,tmaxoffsetindex,shotno)
   
        #pmds.mdsdisconnect()
    #pmds.mdsclose('tcv_shot',shotno)
    connect.closeTree('tcv_shot',shotno)

    axuvstruct.data = (data_raw-offset[i,None])/(0.24*gain*axuvstruct.geometry.etend[:,None])
    
    # 10/2^15:      [Volt/bit] conversion factor of the ADC
    # 0.24:         [A/W] conversion factor according to Boivin's article

    axuvstruct.time=time
    axuvstruct.raw.data=data_raw
    axuvstruct.raw.offset=offset
    axuvstruct.raw.gains=tile(gain,(1,140)).T
    axuvstruct.raw.convfact=0.24
    axuvstruct.raw.shotno_offset=shotno_offset
    axuvstruct.name='axuv'
    axuvstruct.shotno=shotno
    display('AXUV data download has finished.')
    
    return axuvstruct






def axuv_dld(connect, shotno,tbegin=.4,tend=1.8,period=5e-6, diag='PLASMA'):
    # axuv=axuv_dld(shot,tm,tM,period,array);
    # Basic routine to download axuv data
    # Created during 2014 shutdown when crpppc250 was off.
    #
    # Benoit Labit. May 2014.
    

    #geom = axuv_calc_los(axuv_detectorpos('axuv'))
    #names=fieldnames(axuv)
    
    #for i=1:length(names)
        #axuv = setfield(axuv,names{i},[]);
    #end
    #axuv(2)=axuv(1);
    diag_path = os.path.dirname(os.path.realpath(__file__))+'/'
    bolomap=zeros((20,7),dtype=int)
    boards = zeros((20,7),dtype=int)
    lymanmap = zeros((20,7),dtype=int)
    boardsC = zeros((20,7),dtype=int)

    if diag in ('PLASMA', 'BOTH'):
        # Patch panel
        bolomap[:,0]=[1,20,2,19,3,18,4,17,5,16,6,15,7,14,8,13,9,12,10,11]
        bolomap[:,1]=[43,62,44,61,45,60,46,59,47,58,48,57,49,56,50,55,51,54,52,53]
        bolomap[:,2]=[1,94,2,93,3,92,4,91,5,90,6,89,7,88,8,87,9,86,10,85]
        bolomap[:,3]=[43,42,44,41,45,40,46,39,47,38,48,37,49,36,50,35,51,34,52,33]
        bolomap[:,4]=[85,84,86,83,87,82,88,81,89,80,90,79,91,78,92,77,93,76,94,75]
        bolomap[:,5]=[21,42,22,41,23,40,24,39,25,38,26,37,27,36,28,35,29,34,30,33]
        bolomap[:,6]=[65,84,66,83,67,82,68,81,69,80,70,79,71,78,72,77,73,76,74,75]
        
        boards[:,r_[0:2,3:7]]= [11 ,11, 12 ,12 ,13, 13]
        boards[:,2]=[12,11,12,11,12,11,12,11,12,11,12,11,12,11,12,11,12,11,12,11]
        if shotno>=44499:  #the cabling changed after the recommissioning!!!
            bolomap[:,2:5]=flipud(bolomap[:,2:5])
            boards[:,2:5]=flipud(boards[:,2:5])
            axuv1=axuv_calc_los(axuv_detectorpos('axuv'))
            axuv2=axuv_calc_los(axuv_detectorpos('calibration'))

            gain=400e3
            if os.path.isfile(diag_path+'axuv_finiteetend_axuv.mat'):
                geom = loadmat(diag_path+'axuv_finiteetend_axuv.mat') 
                #change is lower than 3%
                axuv1.etend=squeeze(geom['etendstruct']['etend'].item())
        
            
        elif shotno<44499 and shotno>=34666:
            axuv1=axuv_calc_los(axuv_detectorpos('axuv_old'))
            axuv2=axuv_calc_los(axuv_detectorpos('calibration_old'))
            gain=20e3
        elif  shotno>=50000 :#??
            bolomap[:,2] = lymanmap[:,2];
            bolomap[:,3] = lymanmap[:,3];
            bolomap[:,4] = lymanmap[:,4];
            lymanmap[:,2] = bolomap[:,2];
            lymanmap[:,3] = bolomap[:,3];
            lymanmap[:,2] = bolomap[:,2]; #BUG??
            boards[:,2] = boardsC[:,2];
            boards[:,3] = boardsC[:,3];
            boards[:,4] = boardsC[:,5];
            boardsC[:,2] = boards[:,2];
            boardsC[:,3] = boards[:,3];
            boardsC[:,4] = boards[:,4];
            
            axuv1=axuv_calc_los(axuv_detectorpos('axuv'));
            axuv2=axuv_calc_los(axuv_detectorpos('calibration'));
            gain=400e3;
            
        else:
            raise Exception('Too old shot!!')
        
    if diag in ('CALIBRATION', 'BOTH'):
        
        lymanmap[:,0]=[21,42,22,41,23,40,24,39,25,38,26,37,27,36,28,35,29,34,30,33];
        lymanmap[:,1]=[65,84,66,83,67,82,68,81,69,80,70,79,71,78,72,77,73,76,74,75];
        lymanmap[:,2]=[21,20,22,19,23,18,24,17,25,16,26,15,27,14,28,13,29,12,30,11];
        lymanmap[:,3]=[65,62,66,61,67,60,68,59,69,58,70,57,71,56,72,55,73,54,74,53];
        lymanmap[:,4]=[11,10,12,9,13,8,14,7,15,6,16,5,17,4,18,3,19,2,20,1];
        lymanmap[:,5]=[43,62,44,61,45,60,46,59,47,58,48,57,49,56,50,55,51,54,52,53];
        lymanmap[:,5]=[85,64,86,63,87,32,88,31,89,96,90,95,91,64,92,63,93,32,94,31];
        boardsC[:,:6]=[11,11,12,12,13,13]

        boardsC[:,6]=[13,12,13,12,13,12,13,12,13,11,13,11,13,11,13,11,13,11,13,11];

        if shotno>=44499:  #the cabling changed after the recommissioning!!!
            lymanmap[:,2:5]=flipud(lymanmap[:,2:5])
            boardsC[:,2:5]=flipud(boardsC[:,2:5]);
            axuv2=axuv_calc_los(axuv_detectorpos('calibration'))
            gain=400e3;
            
        elif shot<44499 and shot>=34666:
            gain=20e3;
            axuv2=axuv_calc_los(axuv_detectorpos('calibration_old'))
        else:
            raise Exception('The routine is not suitable for downloading these data.');
        
    
    #
    if diag in 'PLASMA':
        BOARD = boards.flatten(order='F')
        CHANNEL = bolomap.flatten(order='F')
    elif diag in 'CALIBRATION':
        BOARD = boardsC.flatten(order='F')
        CHANNEL = lymanmap.flatten(order='F')
    else:
        BOARD=r_[boards.flatten(order='F'), boardsC.flatten(order='F')]
        CHANNEL=r_[bolomap.flatten(order='F'), lymanmap.flatten(order='F')]

    #pmds.mdsopen('tcv_shot',shotno)
    connect.openTree('tcv_shot',shotno)

    nchn=len(CHANNEL)
    #print('%03d/%03d downloaded \r'%(0,nchn))
    tvec=array(connect.get(r'dim_of(\atlas::dt196_axuv_001:CHANNEL_001)'))

    #tvec = (1,1)
    #from mds_data import mds_load
    #nodes = 

    #nodes = 
    
    #mds_load(connect, nodes,chns, chnfmt=None, tmin=-infty,tmax=infty,step=1,raw=True,
             #remove_offset=False,offset_time_min = -infty, offset_time_max=0):
    
    
    
    V_max = 10.
    bitres = 16
    
    convfact = 2*V_max/(2**bitres-2.)
    
    raw = empty((len(tvec),nchn),dtype='int16')
    
    
    #tqdm.trange
    for i in tqdm.trange(nchn):
        
        cms = r'raw_of(\atlas::dt196_axuv_%.3d:CHANNEL_%.3d)'%(BOARD[i]%10,CHANNEL[i])

        raw[:,i]= connect.get( cms)

        #print('%0.3d/%0.3d downloaded \r'%(i,nchn))
    #import IPython
    #IPython.embed()
 
    #tvec=mdsvalue('dim_of('+cms+')')

    #pmds.mdsclose('tcv_shot',shotno)
    connect.closeTree('tcv_shot',shotno)

    # Remove offset
    #BUG why is not used offset from the data?? 
    #offset = loadmat('/home/labit/matlab/axuv/axuv_offset')
    #if  diag == 'PLASMA':
        #raw-= offset['offset']
    #if diag == 'CALIBRATION':
        #raw-= offset['offsetC']
    #if diag == 'BOTH':
        #raw-= r_[offset['offset'],offset['offsetC']]
    
    # Calibration
    etend = axuv1.etend

     
    data = raw/single(0.24*gain*etend/convfact)
    saturated = raw == (2**bitres-2)/2
    if any(tvec< 0):
        data-=  data[tvec<0].mean(0)
    else:
        offset = loadmat(diag_path+'axuv_offset')
        
        if  diag == 'PLASMA':
            data-= offset['offset']/(0.24*gain*etend)
        if diag == 'CALIBRATION':
            data-= offset['offsetC']/(0.24*gain*etend)
        if diag == 'BOTH':
            data-= r_[offset['offset'],offset['offsetC']]/(0.24*gain*etend)
        
        #import IPython
        #IPython.embed()
        
    #geom1.geometry=tmp(1)
    #axuv(2).geometry=tmp(2);

    if diag == 'PLASMA':
        axuv1.data=data#[:,:140] #BUG can be dimension higher?? 
        axuv1.time=tvec
        axuv1.saturated = saturated
        #axuv1.raw=raw#[:,:140]
        return axuv1

    
    if diag == 'CALIBRATION':
        axuv2.data=data#[:,:140]
        axuv2.time=tvec
        axuv2.saturated = saturated

        return axuv2

    
    if diag == 'BOTH':
        axuv1.data=data[:,:140];
        axuv1.time=tvec;
        #axuv1.raw.data=raw[:,:140];
        axuv2.data=data[:,140:280];
        axuv2.time=tvec;
        axuv1.saturated=saturated[:,:140];
        axuv2.saturated=saturated[:,140:280];

        #axuv2.raw.data=x[:,140:280];
        return axuv1, axuv2





class bolostruct_class:
    def __init__(self):
        class raw:pass
        self.raw = raw()
        class special:pass
        self.special = special()
        
            

    




def get_chord_offset(connect,board,channel,tbeginindex,tendindex,tminoffsetindex,
                     tmaxoffsetindex,shotno):
    
    #connect.openTree('crpppc250.epfl.ch::axuv,$',shotno)
    #import IPython
    #IPython.embed()
    if tbeginindex is None and tendindex is None:
        data=asarray(connect.get('.BOARD_%.2d:CHANNEL_%.3d'%(board,channel)))  #load raw data
        if len(data) > 1:
            offset=mean(data[tminoffsetindex:tmaxoffsetindex])
      
    else:
        data = asarary(connect.get('.BOARD_%.2d:CHANNEL_%.3d[%d:%d]'%(board,channel,tbeginindex,tendindex)))   #load raw data
        if len(data) > 1:
            offset=mean(mdsvalue('.BOARD_%.2d:CHANNEL_%.3d[%d:%d]'%(board,channel,tminoffsetindex,tmaxoffsetindex)))
        
        
    if shotno>=45523 and board==11 and channel==48: 
        data[data<-5]+= 10
        if tbeginindex is None:
            if len(data) > 1:
                offset=mean(data[tminoffsetindex:tmaxoffsetindex])
        else:
            if len(data) > 1:
                offset=mean(mdsvalue('.BOARD_%.2d:CHANNEL_%.3d[%d:%d]'%(board,channel,tminoffsetindex,tmaxoffsetindex)))
    
    data   = data  *(10./2**16/2)
    offset = offset*(10./2**16/2)
    
    return  data,offset



#def get_chord(board, channel,tbeginindex,tendindex):
    #if isempty(tbeginindex):
        #b=mdsvalue( '.BOARD_%.2d:CHANNEL_%.3d'%(board,channel))*(10/2**16/2)
    #else:
        #b=mdsvalue('.BOARD_%.2d:CHANNEL_%.3d[%d:%d]'%(board,channel,tbeginindex,tendindex))*(10/2**16/2)
    #return b



#def (board, channel,tminoffsetindex,tmaxoffsetindex):
    #ofget_offsetfset=mean(mdsvalue('.BOARD_%.2d:CHANNEL_%.3d[%d:%d]'%(board,channel,tminoffsetindex,tmaxoffsetindex)))*(10/2**16/2)

    #return offset
    
    

    
class geometry_struct_class:
    def __init__(self):
        pass

    
    
def  axuv_calc_los(geometric_data):
    #Simple etendue and line of sight calculation program


    xdet=geometric_data.xdet
    ydet=geometric_data.ydet
    d1=geometric_data.d1
    d2=geometric_data.d2
    b1=geometric_data.b1
    b2=geometric_data.b2
    xap=geometric_data.xap
    yap=geometric_data.yap
    vangle=geometric_data.vangle
    detangle=geometric_data.detangle
    slitdist=geometric_data.slitdist
    slitwidth=geometric_data.slitwidth

    ymin=min(geometric_data.vessel_y)
    ymax=max(geometric_data.vessel_y)
    xmin=min(geometric_data.vessel_x)
    xmax=max(geometric_data.vessel_x)

    chordangle=arctan2(yap-ydet,xap-xdet)

    #---calculation of approximate etendues---
    r=hypot(ydet-yap, xdet-xap)
    
    etend=d1*slitwidth*b1*b2*cos(chordangle-vangle)*cos((chordangle-vangle)
            -(detangle-vangle))/r/(r-cos(chordangle-vangle)*slitdist)/4/pi # we have then m^2


    #---calculation of lines of sights---
    m=(ydet-yap)/(xdet-xap)
    b=ydet-m*xdet

    nl=len(xdet)
    xchord=zeros((2,nl))
    ychord=zeros((2,nl))

    xchord[0]=xdet
    ychord[0]=ydet

    iup=vangle<0
    isi=vangle==pi
    ido=vangle>0

    if any(iup):
        ychord[1,iup]=ymin
        xchord[1,iup]=(ychord[1,iup]-b[iup])/m[iup]
    
    if any(ido):
        ychord[1,ido]=ymax
        xchord[1,ido]=(ychord[1,ido]-b[ido])/m[ido]
    
    if any(isi):
        xchord[1,isi]=xmin
        ychord[1,isi]=m[isi]*xchord[1,isi]+b[isi]
    

    ileft= xchord[1]<xmin
    xchord[1,ileft]=xmin
    ychord[1,ileft]=m[ileft]*xchord[1,ileft]+b[ileft]

    irig=xchord[1]>xmax
    xchord[1,irig]=xmax
    ychord[1,irig]=m[irig]*xchord[1,irig]+b[irig]

    iuplong = ychord[1]>ymax
    ychord[1,iuplong]=ymax
    xchord[1,iuplong]=(ychord[1,iuplong]-b[iuplong])/m[iuplong]

    idownlong= ychord[1]<ymin
    ychord[1,idownlong]=ymin
    xchord[1,idownlong]=(ychord[1,idownlong]-b[idownlong])/m[idownlong]
    #---calculation of lines of sight---

    #---creating the psi toolbox object---
    #if ~isempty(which('psitbxdcd'))
        #rd = xdet;
        #zd = ydet;
        #phid = 0;
        #pvd = -(chordangle-pi);
        #pvd(find(pvd>pi))=pvd(find(pvd>pi))-2*pi;
        #pvd(find(pvd<=-pi))=pvd(find(pvd<=-pi))-2*pi;
        #tvd = 0;
        #dcd = psitbxdcd(rd,zd,phid,pvd,tvd);
    #else
        #dcd=[];
    #end
    #import IPython
    #IPython.embed()
    
    #from matplotlib.patches import Rectangle
    #rec = Rectangle([geometric_data.vessel_x[0],geometric_data.vessel_y[0]],
              #diff(geometric_data.vessel_x), diff(geometric_data.vessel_y),facecolor='none')
    
    
    #gca().add_patch(rec)
    #plot(xchord,ychord,lw=.1, c='k')
    #show()
    
    geometry_struct = geometry_struct_class()
    geometry_struct.xchord=xchord
    geometry_struct.ychord=ychord
    geometry_struct.etend=etend
    geometry_struct.chordangle=chordangle
    geometry_struct.vessel_x=geometric_data.vessel_x
    geometry_struct.vessel_y=geometric_data.vessel_y
    geometry_struct.majorcamlim=geometric_data.majorcamlim
    geometry_struct.minorcamlim=geometric_data.minorcamlim
    geometry_struct.detectorpos=geometric_data
    geometry_struct.detectortype=geometric_data.detectortype
    #geometry_struct.dcd=dcd;
        
    return geometry_struct
    
    
    
class geometric_data_class:
    def __init__(self):
        pass
    
def axuv_detectorpos(det_type):
    #[geometric_data]=axuv_detectorpos(det_type)
    #det_type possibilities:
    #       'axuv'
    #       'calibration'
    #       'BOLO'
    #       'axuvproto'

    geometric_data = geometric_data_class()
    
    if det_type in ('axuv','calibration'):
            d1a=0.00075                                     # detector poloidal size
                                                                #According to: http://www.ird-inc.com/axuvarr.html
            d2a=0.004                                       # detector toroidal size
                                                                #According to: http://www.ird-inc.com/axuvarr.html
                                                            #Warning! The effective surface is only approximated by these parameters and should be measured)
            det_length = 0.00095                           #distance between diode centers
            vangle_cp=[-pi/2, -pi/2,pi, pi, pi, pi/2, pi/2]       # angle of pinhole surface normal to the horizontal
                                                                    #the direction is: from the detector to the pinhole

            num_channels = 20                                   #number of diodes per detector
            xpos=[0.719, 1.04, 1.20, 1.20, 1.20, 1.04, 0.719]                    # x positions of the pinholes
            ypos=[0.795, 0.795, 0.4526, -0.0025, -0.4576, -0.795, -0.795]        # y positions of the pinholes
                    #---According to A. Degeling...---
            geometric_data.majorcamlim= r_[20.5:140:20]                                     
            geometric_data.minorcamlim=[]
            ae_design=[-0.004,0.005,-0.0024,0,0.0024,-0.005,0.004]         #parallel shift of the detector array relative to the pinhole
            da_design=[0.012,0.015,0.009,0.0078,0.009,0.015,0.012]         #pinhole-array perpendicular poloidal distance


            slitwidtha=1e-3                                                        #slit toroidal size
            slitheighta=1e-3                                               #slit perpendicular size
            if det_type == 'axuv':      #this is the left camera
                b1c=array([0.48,0.49,0.49,0.465,0.49,0.48,0.51])*1e-3             # aperture poloidal size            
                b2c=array([0.97,0.97,0.98,0.97,0.99,0.97,1])*1e-3                 # aperture toroidal size

                #load axuv_detectorpos_laser            
                #load axuv_detectorpos_gauge
                #slitdistc=sens_slit_dist[0]

                slitdistc     =array([.6125,.6825,.6500,.7275,.6975,.7625,.8100])*1e-3    #according to gauge measurements
                del_a_par_left=array([.0695,.1244,.0222,-.2831,-.0510,.1135,.1954])*1e-3  #according to laser measurements
                del_d_par_left=array([.1741,.2223,.2234,.3190,.2381,.2326,.1320])*1e-3     #according to laser measurements
                ae=ae_design+del_a_par_left
                da=da_design+del_d_par_left
            elif  det_type == 'calibration':   #this is the right camera
                b1c=array([0.495,0.47,0.5,0.47,0.5,0.485,0.55])*1e-3              # aperture poloidal size            
                b2c=array([0.99,0.97,0.98,0.97,0.99,0.98,1.1])*1e-3               # aperture toroidal size
                #load axuv_detectorpos_laser
                #load axuv_detectorpos_gauge
                #slitdistc=sens_slit_dist(2,:);     #according to measurement with the gauge
                slitdistc=array([.5850,.7125,.6750,.6925,.6050,.6750,.8700])*1e-3                  #according to gauge measurements
                del_a_par_right=array([.0399,.2102,-.1040,-.2470,-.2460,-.0599,.2010])*1e-3        #according to laser measurements
                del_d_par_right=array([.1935,.2547,.2698,.2444,.2162,.2254,.2074])*1e-3            #according to laser measurements
                ae=ae_design+del_a_par_right
                da=da_design+del_d_par_right
            
    elif det_type in ('axuv_old','calibration_old'):
        d1a=0.00075                                     # detector poloidal size
                                                            #According to: http://www.ird-inc.com/axuvarr.html
        d2a=0.004                                       # detector toroidal size
                                                            #According to: http://www.ird-inc.com/axuvarr.html
                                                        #Warning! The effective surface is only approximated by these parameters and should be measured)
        det_length = 0.00095                           #distance between diode centers
        vangle_cp=[-pi/2,-pi/2,pi,pi,pi,pi/2,pi/2]       # angle of pinhole surface normal to the horizontal
                                                                #the direction is: from the detector to the pinhole

        num_channels = 20                             #number of diodes per detector
        xpos=[0.719,1.04,1.20,1.20,1.20,1.04,0.719]                    # x positions of the pinholes
        ypos=[0.795,0.795,0.4526,-0.0025,-0.4576,-0.795,-0.795]        # y positions of the pinholes
                #---According to A. Degeling...---
        geometric_data.majorcamlim=r_[20.5:140,20]                                       
        geometric_data.minorcamlim=[]
        ae_design=[-0.004,0.005,-0.0024,0,0.0024,-0.005,0.004]         #parallel shift of the detector array relative to the pinhole
        da_design=[0.012,0.015,0.009,0.0078,0.009,0.015,0.012]         #pinhole-array perpendicular poloidal distance

        b1c=ones(7)*0.5e-3
        b2c=ones(7)*5e-3

        slitdistc=zeros_like(b2c)
        slitwidtha=d2a                                 #slit toroidal size
        slitheighta=0                                  #slit perpendicular size

        ae=ae_design
        da=da_design
        
    elif det_type == 'axuvproto':
        raise Exception('Not implemented yet.')
        d1a=0.002;                                       # detector poloidal size
                                                            #According to: http://www.ird-inc.com/axuvarr.html
        d2a=0.005;                                       # detector toroidal size
                                                        #According to: http://www.ird-inc.com/axuvarr.html
                                                        #Warning! The effective surface is only approximated by these parameters and should be measured)
        det_length = 0.00212;                             #distance between diode centers

        vangle_cp=[-pi/2,pi/2];                            # angle of pinhole surface normal to the horizontal
                                                                #the direction is: from the detector to the pinhole
        num_channels = 16;                              #number of diodes per detector
        #---According to A. Degeling: The AXUV Bolometer and La Camera
        #System on TCV, 4th Chapter: Storage of the camera parameters
        #filename: xb_c3_d16u2.mat and xb_c3_d16u2.m
        #original place: hal.epfl.ch:22/home/degeling/bolometry
        #detector type: AXUV16ELO
        #gain = [1.2,1.2].*1e6          #Pre amplifier gains (ohms) (what does it mean????)
        xpos=[0.831,1.0395];          # x position of the pinholes
        ypos=[0.916,-0.8185];          # y position of the pinholes
        b1c=[0.0005,0.00056];                  # aperture poloidal size
        b2c=[0.005,0.005];                      # aperture toroidal size

        ae=[0.007,-0.005];                   #parallel shift of the detector array relative to the pinhole
        da=[0.07,0.0959];                  # pinhole-array perpendicular poloidal distance
        #---According to A. Degeling...---
        geometric_data.majorcamlim=[16.5]
        geometric_data.minorcamlim=[]

    elif det_type in  'BOLO':
        vangle_cp=[-pi/2,pi,pi,pi,pi,pi,pi,pi/2]                            # angle of pinhole surface normal
                                                                        #zero at nine o'clock 
                                                                        #direction is inside the chamber
                                                                        #surprisingly counterclockwiseclockwise 

        xpos=[0.88,1.235,1.235,1.235,1.235,1.235,1.235,0.88]       # x position of the pinholes
        ypos=[0.815,0.455,0.455,-.0025,-.0025,-0.46,-0.46,-0.815]          # y position of the pinholes
        #Bernard Joye's original data: geometric_data.mat

        # position of the detectors
        xdet=[0.862,0.867,0.8721,0.8771,0.8829,0.8879,0.893,0.898,
            1.278,1.2806,1.2833,1.2861,1.2883,1.2893,1.2902,1.2912]
        xdet[16:24]=xdet[8:16][::-1]
        xdet[24:40]=xdet[8:24]
        xdet[40:56]=xdet[8:24]
        xdet[56:64]=xdet[0:8][::-1]

        ydet=[0.9057,0.9063,0.9069,0.9075]
        ydet[4:8]=ydet[:4][::-1]
        ydet[8:24]=[0.4916,0.4873,0.483,0.4787,0.473,0.468,0.4631,0.4581,\
            0.4519,0.4469,0.442,0.437,0.4313,0.427,0.4227,0.4184]
        ydet[24:40]=array(ydet[8:24])-0.4575
        ydet[49:56]=array(ydet[8:24])-0.915
        ydet[56:64]=-array(ydet[0:8])
        #xdet and ydet are checked, they are the same but with better rounding as
            #\base::bolo:radial_pos and \base::bolo:z_pos 
        
        
        d1a=0.0015                                  # detector poloidal size
        d2a=0.004                                    # detector toroidal size
        b1c=array([0.0026,0.0022,0.0022,0.0022,0.0022,0.0022,0.0022,0.0026])*2           # aperture poloidal size

        b2c=array([0.01,0.008,0.008,0.008,0.008,0.008,0.008,0.01])*2                     # aperture toroidal size
        slitdistc=zeros_like(b2c)
        slitwidtha=d2a 
        slitheighta=0

        num_channels= 8                             #number of detector channels in one camera
        num_channelplane=4                            #number channels in one detector array

        geometric_data.majorcamlim=[8.5,24.5,40.5,56.5]
        geometric_data.minorcamlim=[4.5,12.5,16.5,20.5,28.5,32.5,36.5,44.5,48.5,52.5,60.5]


    #---calculating some properties of each chords---
    geometric_data.b1=repeat(b1c, num_channels)
    geometric_data.b2=repeat(b2c, num_channels)

    geometric_data.xap=repeat(xpos, num_channels)
    geometric_data.yap=repeat(ypos, num_channels)

    geometric_data.vangle= repeat(vangle_cp, num_channels)

    geometric_data.d1=ones_like(geometric_data.b1)*d1a
    geometric_data.d2=ones_like(geometric_data.b1)*d2a

    #slitdist=ones(num_channels,1)*slitdistc
    geometric_data.slitdist=repeat(slitdistc, num_channels)
    geometric_data.slitwidth= repeat(slitwidtha,len(geometric_data.vangle))
    
    geometric_data.slitheight= repeat(slitheighta, num_channels)
    #---calculating some properties of each chords---

    #---calculation of detector positions and detector angles---
    if  det_type in ('axuv','calibration','axuv_old','calibration_old','axuvproto'):
        #---neccessary because the the relative position of the diodes to the
        #   aperture is given---
        aee = repeat(ae,num_channels)
        #print shape(da)
        rdio1=repeat(-array(da)*cos(vangle_cp), num_channels)
        rdio1= rdio1 +geometric_data.xap
        zdio1= repeat(-array(da)*sin(vangle_cp), num_channels)
        zdio1= zdio1+geometric_data.yap 

        indd=arange(num_channels)-(num_channels-1)/2.        #indices of diodes in one detector
        delta_dio=indd*det_length
        ddelta_dio= tile(delta_dio, len(vangle_cp))
    
        ddelta_dio=ddelta_dio+aee
        geometric_data.xdet=rdio1-ddelta_dio*sin(geometric_data.vangle)       #x <r> coordinates of diodes
        geometric_data.ydet=zdio1+ddelta_dio*cos(geometric_data.vangle)       #y <z> coordinates of diodes
        #---neccessary because the the relative position of the diodes to the
        #   aperture is given---
        geometric_data.detangle=geometric_data.vangle
        #geometric_data.detangle(find(geometric_data.detangle>=pi))=geometric_data.detangle(find(geometric_data.detangle>=pi))-2*pi;
    elif det_type in 'BOLO':
        geometric_data.detangle = zeros_like(geometric_data.vangle)
        for i in arange(0, len(geometric_data.vangle),num_channelplane):
            geometric_data.detangle[i:i+num_channelplane]=arctan2(-(xdet[i+num_channelplane-1]-xdet[i]),ydet[i+num_channelplane-1]-ydet[i])#-vangle(i:i+num_channelplane-1);
        
        geometric_data.xdet=xdet
        geometric_data.ydet=ydet
    
    #---calculation of detector positions and detctor angles---

    geometric_data.vessel_y=[-0.75, 0.75]
    geometric_data.vessel_x=[0.624, 1.136]
    geometric_data.detectortype=det_type


    return geometric_data
    
    
    
    
    

#geometric_data = axuv_detectorpos('axuv')
#axuv_calc_los(geometric_data)
#geometric_data = axuv_detectorpos('calibration')
#axuv_calc_los(geometric_data)
#geometric_data = axuv_detectorpos('BOLO')
#axuv_calc_los(geometric_data)

#axuv_detectorpos('axuvproto')



#exit()
#det_type possibilities:
#       'axuv'
#       'calibration'
#       'BOLO'
#       'axuvproto'


#'gottardi', 'bessel', 'golay', 'interpos', 'butter'
#server =  'localhost:8002'
#pmds.mdsconnect( server)
##shotno = 40412
#axuv_dld(45186)
##get_bolo(shotno,filtertype='golay')

#exit()

    
#server =  'localhost:8002'
#pmds.mdsconnect( server)
#shotno = 40412
#axuv_get_axuv(shotno)
    
    
    
    #server =  'localhost:8003'
#pmds.mdsconnect( server)
    
#mdsconnect( 'tcvdata.epfl.ch')  
#[axuvstruct]=axuv_get_axuv(45104,0.2,0.3);
#[axuvstruct]=axuv_get_axuv(45190,0.2,0.3);
#[axuvstruct]=axuv_get_axuv(45186,0.2,0.3);

#mdsdisconnect


