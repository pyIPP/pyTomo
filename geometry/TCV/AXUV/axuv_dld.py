
from numpy import *
#import pmds
from matplotlib.pylab import *
#import ././tqdm
import os.path 

import sys,os
#print os.path.abspath('../../')
sys.path.append(os.path.abspath('../../'))
#import tqdm
from scipy.io import loadmat


from get_axuv import axuv_calc_los, axuv_detectorpos


def axuv_dld(shotno,tbegin=0.4,tend=1.8,period=5e-6,array='PLASMA'):
    # axuv=axuv_dld(shot,tm,tM,period,array);
    # Basic routine to download axuv data
    # Created during 2014 shutdown when crpppc250 was off.
    #
    # Benoit Labit. May 2014.

    axuv =axuv_calc_los(axuv_detectorpos('axuv'));
    #names=fieldnames(axuv);
    #for i=1:length(names)
        #axuv = setfield(axuv,names{i},[]);
    
    #axuv(2)=axuv(1);

    bolomap=[]; boards=[]; lymanmap=[]; boardsC=[];


    # Patch panel "PLASMA" array
    bolomap(:,1)=[1 20 2 19 3 18 4 17 5 16 6 15 7 14 8 13 9 12 10 11];
    bolomap(:,2)=[43 62 44 61 45 60 46 59 47 58 48 57 49 56 50 55 51 54 52 53];
    bolomap(:,3)=[1 94 2 93 3 92 4 91 5 90 6 89 7 88 8 87 9 86 10 85];
    bolomap(:,4)= [43 42 44 41 45 40 46 39 47 38 48 37 49 36 50 35 51 34 52 33];
    bolomap(:,5)=[85 84 86 83 87 82 88 81 89 80 90 79 91 78 92 77 93 76 94 75];
    bolomap(:,6)=[21 42 22 41 23 40 24 39 25 38 26 37 27 36 28 35 29 34 30 33];
    bolomap(:,7)=[65 84 66 83 67 82 68 81 69 80 70 79 71 78 72 77 73 76 74 75];

    boards(:,[1:2 4:7])=repmat([11 11  12 13 13 13],20,1);
    boards(:,3)=[12 11 12 11 12 11 12 11 12 11 12 11 12 11 12 11 12 11 12 11];

    # Patch panel "CALIBRATION" array
    lymanmap(:,1)=[21 42 22 41 23 40 24 39 25 38 26 37 27 36 28 35 29 34 30 33];
    lymanmap(:,2)=[65 84 66 83 67 82 68 81 69 80 70 79 71 78 72 77 73 76 74 75];
    lymanmap(:,3)=[21 20 22 19 23 18 24 17 25 16 26 15 27 14 28 13 29 12 30 11];
    lymanmap(:,4)=[65 62 66 61 67 60 68 59 69 58 70 57 71 56 72 55 73 54 74 53];
    lymanmap(:,5)=[11 10 12 9 13 8 14 7 15 6 16 5 17 4 18 3 19 2 20 1];
    lymanmap(:,6)=[43 62 44 61 45 60 46 59 47 58 48 57 49 56 50 55 51 54 52 53];
    lymanmap(:,7)=[85 64 86 63 87 32 88 31 89 96 90 95 91 64 92 63 93 32 94 31];

    boardsC(:,[1:6])=repmat([11 11 12 12 13 13],20,1);
    boardsC(:,7) =  [13 12 13 12 13 12 13 12 13 11 13 11 13 11 13 11 13 11 13 11];


    if shotno>=44499 && shotno<50000  #the cabling changed after the recommissioning!!!
        bolomap(:,[3 4 5])=flipud(bolomap(:,[3 4 5]));
        boards(:,[3 4 5])=flipud(boards(:,[3 4 5]));
        lymanmap(:,[3 4 5])=flipud(lymanmap(:,[3 4 5]));
        boardsC(:,[3 4 5])=flipud(boardsC(:,[3 4 5]));
        axuv(1)=axuv_calc_los(axuv_detectorpos('axuv'));
        axuv(2)=axuv_calc_los(axuv_detectorpos('calibration'));
        gain=400e3;
    elseif shotno<44499 && shotno>=34666
        axuv(1)=axuv_calc_los(axuv_detectorpos('axuv_old'));
        axuv(2)=axuv_calc_los(axuv_detectorpos('calibration_old'));
        gain=20e3;
    elseif shotno>=50000 
        bolomap(:,3) = lymanmap(:,3);
        bolomap(:,4) = lymanmap(:,4);
        bolomap(:,5) = lymanmap(:,5);
        lymanmap(:,3) = bolomap(:,3);
        lymanmap(:,4) = bolomap(:,4);
        lymanmap(:,3) = bolomap(:,3);
        boards(:,3) = boardsC(:,3);
        boards(:,4) = boardsC(:,4);
        boards(:,5) = boardsC(:,5);
        boardsC(:,3) = boards(:,3);
        boardsC(:,4) = boards(:,4);
        boardsC(:,5) = boards(:,5);
        
        axuv(1)=axuv_calc_los(axuv_detectorpos('axuv'));
        axuv(2)=axuv_calc_los(axuv_detectorpos('calibration'));
        gain=400e3;
        
    else
        error('Too old shot!!')
        
    end

    #
    clear x
    BOARD=[boards(:);boardsC(:)];
    CHANNEL=[bolomap(:);lymanmap(:)];

    # Keep only what you need
    if strcmpi(array,'PLASMA')
        BOARD=BOARD(1:140); CHANNEL=CHANNEL(1:140);
    end
    if strcmpi(array,'CALIBRATION')
        BOARD=BOARD(141:280); CHANNEL=CHANNEL(141:280);
    end


    clear functions
    mdsopen(shotno)
    nchn=length(CHANNEL);
    fprintf('#03d/#03d downloaded \r',0,nchn);

    for i=1:nchn
        tmp=num2str(BOARD(i));tmp=tmp(2);
        str=['\atlas::dt196_axuv_00',tmp,':CHANNEL_0' num2str(CHANNEL(i),'#0.2d'),...
            '[',num2str(tbegin),':',num2str(tend),':',num2str(period),']'];
        x(:,i)=mdsvalue(str);
        str_off=['\atlas::dt196_axuv_00',tmp,':CHANNEL_0' num2str(CHANNEL(i),'#0.2d'),...
            '[',num2str(-0.04),':',num2str(-0.01),':',num2str(1e-4),']'];
        xoffset(:,i)=mdsvalue(str_off);
        fprintf('#03d/#03d downloaded \r',i,nchn);
        
    end
    fprintf('\n')
    t=mdsvalue(['dim_of(',str,')']);

    mdsclose;
    #keyboard
    # Remove offset
    load('/home/labit/matlab/axuv/axuv_offset')
    #x = x -repmat(mean(xoffset),[size(x,1),1]);
    if shotno >=50000
        offset=offset2015;offsetC=offsetC2015;
    end
    if strcmpi(array,'PLASMA')
    x=x-repmat(offset,[size(x,1),1]);
    end
    if strcmpi(array,'CALIBRATION')
    x=x-repmat(offsetC,[size(x,1),1]);
    end
    if strcmpi(array,'BOTH')
    x=x-repmat([offset,offsetC],[size(x,1),1]);
    end
    # Calibration
    etend=cat(2,axuv(:).etend);


    for i=1:nchn
        
        xx(:,i)=x(:,i)/0.24/gain/etend(i);
    #xx(:,i) = x(:,i);
    end
    tmp=axuv; clear axuv;
    axuv(1).geometry=tmp(1);axuv(2).geometry=tmp(2);

    if strcmpi(array,'PLASMA')
        
        axuv(1).data=xx(:,1:140);
        axuv(1).time=t;
        axuv(1).raw.data=x(:,1:140);
        axuv(2).data=[];
        axuv(2).time=[];
        axuv(2).raw.data=[];
        axuv=axuv(1);
    end
    if strcmpi(array,'CALIBRATION')
        axuv(1).data=[];
        axuv(1).time=[];
        axuv(1).raw.data=[];
        axuv(2).data=xx(:,1:140);
        axuv(2).time=t;
        axuv(2).raw.data=x(:,1:140);
        axuv=axuv(2);
    end
    if strcmpi(array,'BOTH')
        axuv(1).data=xx(:,1:140);
        axuv(1).time=t;
        axuv(1).raw.data=x(:,1:140);
        axuv(2).data=xx(:,141:280);
        axuv(2).time=t;
        axuv(2).raw.data=x(:,141:280);
    end

    return axuv
