#!/usr/bin/env python
# -*- coding: utf-8 -*-


#ssh todstrci@lac911.epfl.ch -L 8002:tcv1.epfl.ch:8000
#ssh todstrci@lac7.epfl.ch -L 8002:tcv1.epfl.ch:8000
#ssh todstrci@lac7.epfl.ch -L 8002:tcvdata.epfl.ch:8000



from numpy import *
#import pmds
from matplotlib.pylab import *
#import ././tqdm
 
import sys,os
#print os.path.abspath('../../')
parent_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
#print parent_path
sys.path.append(parent_path)
import tqdm

def mds_load(connect, nodes,chns, chnfmt=None, tmin=-infty,tmax=infty,step=1,raw=True,
             remove_offset=False,offset_time_min = -infty, offset_time_max=0):

    #TODO načítat celé listy nodů?? 

    #[y,t,f_Nyq]=mdsdownload(kvl) - Download data through open mds connection
    #
    #   node .......
    #      '\dt200::dt200_001'           for   dt200_001 [default]
    #      '\hextip::dt196_001'   for   HEXTIP
    #      '\trcf_torpex_00?'     for   slow CAMACS
    #      '\tr3412_torpex_00?'   for   fast CAMACS
    #
    #   chnfmt ..... sprintf format string for channel node names. May be
    #                'auto' for an automatic determination.
    #   chn ........ Vector of channel numbers to download. The time base is
    #                taken from the first element.
    #   tm ......... Minimum time.
    #   tM ......... Maximum time.
    #   inc ........ Step size (e.g. inc = 10: download every 10th point).
    #   raw ........ Download data as 'raw' (int16) and convert locally.
    #                Note that the 'rmo' offset removal will then be performed
    #                using integer arithmetics, meaning the result will be
    #                accurate to 1 bit.
    #   DTYPE ...... 'single' or 'double'
    #   t_units .... Unit of time basis {'s_Hz','ms_kHz','mus_MHz','ns_GHz'}
    #   rmo ........ Remove offset by performing an average over the
    #                after-shot phase before downloading the data.
    #   rmo_offs ... Number of points to wait after the magnetron shutdown
    #                before the after-shot phase starts.
    #   rmo_Lmin ... Minimum number of points in the after-shot phase that have
    #                to be available to perform an offset removal. If there are
    #                less, the "mdsdownload:rmo:TooFewPoints" error is thrown.
    #   rmo_inc .... Sampling increment for the offset (don't use values
    #                corresponding to likely pickup frequencies).
    #
    #   Errors thrown:
    #    - mdsdownload:PackNotFound
    #    - mdsdownload:NoMDSConn
    #    - mdsdownload:NoTimeBase
    #    - mdsdownload:rmo:TooFewPoints
    #
    #   S. H. Mueller, 2005/02/19
    #
    #  More or less updated for TCV by B. Labit (march 2009)
    #
    # Fixed a bugg for saturated signals. B. Labit. April 2010

    #if not chns is list and  not chns is tuple:
        #chns = [chns, ]
    #if not nodes is list and not nodes is tuple:
        #nodes = [nodes, ]

    if chnfmt is None:
        chnfmt = 'channel_%.3d'
        
    
    try:
        #shot = pmds.mdsvalue( '$shot')
        shot = int(connect.get( '$shot'))
        if not isscalar(shot):
            raise Exception('No MDS connection open')
    except:
        #pmds.mdsvalue(
        raise Exception('MDS failured')
    
    #if chn == 'auto':
    #siz = size(chn)

    fmtstr = (nodes[0]+ ':'+ chnfmt).replace('/','//')

    # set up tdi commands to obtain the information necessary to reconstruct
    # the time base

    get_tb = '_t=dim_of('+fmtstr%chns[0][0]+') '

    get_ind = '_a=dim_of(data(_t)); '

    set_ret = '[_a[0],_a[size(_a)-1],slope_of(axis_of(_t)),_t[0]]';

    mdscmd = get_tb+get_ind+set_ret



    # fast way to obtain the time base - compute time vector as
    # t = (t_range(4)+[t_range(1):t_range(2)]*t_range(3))*1e-9;
    #print mdscmd
    #exit()
    #import IPython
    #IPython.embed()
    #_t=dim_of(\ATLAS::dt196_xtomo_001:channel_001); _a=dim_of(data(_t)); [_a[0],_a[size(_a)-1],slope_of(axis_of(_t)),_t[0]]
    #mdscmd = r'_t=dim_of(\ATLAS::dt196_xtomo_001:channel_001); _a=dim_of(data(_t)); [_a[0],_a[size(_a)-1],slope_of(axis_of(_t)),_t[0]]'
    #t_range = pmds.mdsvalue(r'_t=dim_of(\ATLAS::dt196_xtomo_001:channel_001)[0] ; _a=dim_of(raw_of(_t))')
    

    #if t_range is None:
        #raise Exception('Cannot determine time base for shot #%d, node %s',shot,chn[0])
    
    
    #tvec = t_range[3]+arange(t_range[0],t_range[1])*t_range[2]
    #print get_tb
    #tvec = pmds.mdsvalue(get_tb)
    tvec = asarray(connect.get(get_tb))
    #print tvec
    
    
    t_ind = slice(tvec.searchsorted(tmin),tvec.searchsorted(tmax)-1, step)
    #print len(tvec), tvec[0], tvec[t_ind.stop]-tmax
    


    # range string
    if isinf(tmin) and isinf(tmax):
        # if all data is required, it's a little faster to do without
        # range indexation
        irange = ''
    else:
        irange = '[%d:%d:%d]'%(t_ind.start, t_ind.stop-1, t_ind.step)
    #print tmin, tmax, irange 
    #exit()
    #print
    
    if remove_offset:
        
        o_ind = slice(tvec.searchsorted(offset_time_min),
                        tvec.searchsorted(offset_time_max))

        ioffset = '[%d:%d]'%(o_ind.start, o_ind.stop)
        


    #download raw data (is integers)
    #BUG range and bitrange are set pernamently!!!
    V_max = 10.
    bitres = 16
    
    convfact = 2*V_max/(2**bitres-2.)
    
    #raw_of(\ATLAS::dt196_lang_001:channel_001)
    tvec = tvec[t_ind]
    #data = []
    
    from tqdm import tqdm, trange
    #print chn
    n_sigs = [len(chn) for chn in chns]
    Nch = sum(n_sigs)
    data = empty((len(tvec),Nch))

    for ich in trange(Nch):
        i_node = where( cumsum(n_sigs) >ich)[0][0]
        node = nodes[i_node]
        chn = chns[i_node][int(ich-sum(n_sigs[:i_node]))]


        datacmd = '_x='
        datacmd +=  'raw_of' if  raw else 'data'
        fmtstr = (node+ ':'+ chnfmt).replace('/','//')
        datacmd += '('+fmtstr%chn+')'+irange
        
   
        data[:,ich] = asarray(connect.get(datacmd))
        #print sig.shape

        if remove_offset:
            data[:,ich] -= connect.get('word(mean(float(_x%s)))'%(ioffset))
        if data.dtype in (int16, int32):
            data[:,ich] = data[:,ich]*float32(convfact)
        #data[:,ich] = sig
        
    #import IPython
    #IPython.embed()
    lens = array([len(d) for d in data])
    if any(lens != median(lens)):
        print(median(lens), lens)
        raise Exception('not all DAS has provided the same length if signals')
    
    
        
    return tvec,vstack(data)


#server =  'localhost:8002'
#pmds.mdsconnect( server)
#pmds.mdsopen('tcv_shot',40412)
#node = '\ATLAS::dt196_xtomo_002'
#chnfmtm = 'channel_%.3d'

#command = r'_x=raw_of(\ATLAS::dt196_xtomo_001:channel_001)'
#command = r'_x=data(\ATLAS::dt196_xtomo_001:channel_001)'

#command = r'_x=raw_of(\ATLAS::dt196_lang_001:channel_001);_x-float(word(mean(float(_x[0:100]))))'
#command = r'_x=raw_of(\ATLAS::dt196_lang_001:channel_001)'

#from time import time

#data = pmds.mdsvalue( command)
#chn  = range(1,10)
#(r'\base::xtomo', 'array_%.3d', range(1,11)),

        #(r'\base::trcf_bolo_002', 'channel_%.3d', range(1,17)),
        #(r'\base::trcf_bolo_003', 'channel_%.3d', range(1,17)),
        #(r'\base::trcf_bolo_004', 'channel_%.3d', range(1,17)),]



#for node, chform, chn in bolo:
    #tvec, data = mds_load(node,chn, chform,remove_offset=True)
    #print tvec.shape, data.shape
    #plot( tvec, data);show()

#t  = time()
#xtomo2 = ( (r'\base::xtomo', 'array_%.3d', range(1,10)),)
#for node, chform, chn in xtomo2:
    #tvec, data = mds_load(node,chn, chform,remove_offset=True)
    #print tvec.shape, data.shape
    ##plot( data.mean(0));show()
#print time()-t

#t  = time()

#for node, chform, chn in xtomo:
    #tvec, data = mds_load(node,chn, chform,remove_offset=True,raw=False)
    #print tvec.shape, data.shape
    ##plot( data.mean(0));show()
#print time()-t
#t  = time()



    
#pmds.mdsdisconnect()












##exit()


###bolo = [(r'\base::trcf_bolo_%.3d'%i, 'channel_%.3d', range(1,17)) for i in range(1,5)]


#data = [(r'\atlas::dt100_northeast_001', 'channel_%.3d',range(1,33)),
        #(r'\atlas::dt100_northeast_002', 'channel_%.3d',range(1,33)),
        #(r'\atlas::dt100_northeast_003', 'channel_%.3d',range(1,33)) ]

'\\ATLAS::DT100_NORTHEAST_001:CHANNEL_01[*:*:*]'

#data = [
 #(r'\atlas::dt196_axuv_011', 'channel_%.3d', range(1,97)),
 #(r'\atlas::dt196_axuv_012', 'channel_%.3d', range(1,97)),]


#data = [
 #(r'\atlas::dt196_xtomo_001', 'channel_%.3d', range(1,97)),
 #(r'\atlas::dt196_xtomo_002', 'channel_%.3d', range(1,97)),
 #(r'\atlas::dt196_xtomo_003', 'channel_%.3d', range(1,9))]



#DMPX = [(r'\atlas::dt100_northeast_001:selected', 'channel_%.3d',range(1,33)),
        #(r'\atlas::dt100_northeast_002:selected', 'channel_%.3d',range(1,33)),
        #(r'\atlas::dt100_northeast_003:selected', 'channel_%.3d',range(1,33)) ]

#server =  'localhost:8003'
###server = 'tcvdata.epfl.ch:8000'
###server = 'tcv1.epfl.ch'

#pmds.mdsconnect( server)
#pmds.mdsopen('tcv_shot',40412)

#for node, chform, chn in data:
    #tvec, data = mds_load(node,chn, chform,remove_offset=True)

    #plot(data.mean(0))
#show()

####tstart=mdsvalue(r'\atlas::dt196_xtomo_001:trig_in') #trigger
####dt  = mdsvalue('\atlas::dt196_xtomo_001:clock_period')*1e-6  #time resolution


#xtomo = [
 #(r'\atlas::dt196_xtomo_001', 'channel_%.3d', range(1,97)),
 #(r'\atlas::dt196_xtomo_002', 'channel_%.3d', range(1,97)),
 #(r'\atlas::dt196_xtomo_003', 'channel_%.3d', range(1,9))]


#for node, chform, chn in xtomo:
    #tvec, data = mds_load(node,chn, chform,remove_offset=True)
    #print tvec.shape, data.shape
    ##plot( data.mean(0));show()
####print time()-t

#exit()

#server =  'tcvdata.epfl.ch:8002'
#pmds.mdsconnect( server)
#xtomo = [
 #(r'\atlas::dt196_axuv_011', 'channel_%.3d', range(1,97)),
 #(r'\atlas::dt196_axuv_012', 'channel_%.3d', range(1,97)),]
#pmds.mdsopen('tcv_shot',45186)

#for node, chform, chn in xtomo:
    #tvec, data = mds_load(node,chn, chform,remove_offset=True)













