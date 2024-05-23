#from loader import * 
import os
import MDSplus as mds
import time
#TODO lod only a data slice
#TODO calibration for all diagnostics 
#TODO geometry for all diags
from collections import OrderedDict
try:
    import config
except:
    import  config_loc  as config

    pass
from numpy import *
from matplotlib.pylab import *



#def check(shot):
    ##fastest check if the shotfile exist
    ##BUG 
    #status = True

    #return status

verbose = False



#def mds_load((mds_server,tree, shot,  TDI)):
    #MDSconn = mds.Connection(mds_server )
    #MDSconn.openTree(tree, shot)
    #data = [MDSconn.get(tdi).data()f for tdi in TDI]
    #MDSconn.closeTree(tree, shot)

    ##print mds_server,tree, shot,  TDI

    #return data




#Rvert(1:7) = 230.4;
#Zvert(1:7) = 73.1;
#Rvert(8:12) = 227.2;
#Zvert(8:12) = 78.1;
#Rvert = Rvert/100;  % convert all units to [m]
#Zvert = Zvert/100;



#;SXR Toroidal Array Geometry for 2009
#;
	#rs 	= 2.3403	;real: 2.3361 grid r of slit on gato grid of sxr camera: avg of three slits
	#zs 	= 0.7856	;real: 0.7862 grid z of slit on sxr camera
	#aa	= 12.0E-6	;slit area
	#ha  = 1.0E-3		;height of slit
	#wa  = 12.0E-3		;width of slit
	#;DETECTOR parameters (AXUV-20ELG from IRD) (MKS units)
	#nae   	= 20		;number of active elements (<=12)
	#nc	= 12		;number of channels
	#vc	= [1,0,1,0,1,0,1,0,1,1,1,1,1,1,1,0,1,0,0,0] ;array identifying the valid channels
  #;vc	= [0,0,0,0,1,0,1,0,1,1,1,1,1,1,1,0,1,0,0,0] ;array identifying the valid channels  remove 1 and 2
	#ec	= [1,1,1,1,1,1,1,1,1,1,1,1]			;optimize for p against these channels
	#rd1 = 2.3364 - 0.00	;r of detector end  	
	#zd1 = 0.8182 + 0.002	;z of detector end 
	#rd2 = 2.3676 - 0.0005		
	#zd2 = 0.8016 + 0.0055	
	#ae  = 3.0E-6		;active area (m^2)
	#he  = 1.0E-3		;element height
	#we  = 3.0E-3		;element width
	#dbd = 5.8964E-3 	;distance from rd1,zd1 to first element
	#db  = 1.05E-3  		;distance between elements
	#dso=2.54E-2		;measured distance between detector and slit
	#nstep = 240		;number of divisions of each path
	#nsubc = 100		;number of subchords for Gaussian weighting scheme
	#nsub = 600		;number of points along subchord
	#nsubl = 513
	#resp = 0.26*1.0E6 	;detector responsisivity uA/W
	#rterm = 0.5		;terminating resistor divides signal by half




def xraypaths(shot, toroidal=False):

    if toroidal:

        #------------------------------------------------------------------------------
        # Toroidal SXR for shot >= 91000
        # SX45R1S1-12 (Toroidal X-ray), DIII-D shot > 91000 (since '97)
        # 98/05/11 from Snider. Took from EFIT /link/efit/sxrt97.dat
        # 98/05/19 modified by Snider.
        # 2002/01/21 added new sxr locations for shot >= 108690 per G.Jackson
        #------------------------------------------------------------------------------
        if shot >= 177100:
            #from 2019  
            #sxrgeom_plot3: load in geometry for 2017 SXR45U view chords
            zslit_45= array(( 73.1,)*8+( 78.1,)*4)/100
            rslit_45= array((230.4,)*8+(227.2,)*4)/100
            
            rwall_45 = array([2.19 , 1.98 ,  1.53, 1.12  , 1.01 , 1.01 , 1.01 ,1.01, 1.01, 1.01, 1.01, 1.01]);
            zwall_45 = array([-0.83 ,-1.06 ,-1.24 ,-1.34 ,-1.01 ,-0.65 ,-0.28 ,0.06, 0.20, 0.37, 0.61, 0.81]);

            xangle_45 = rad2deg(arctan2(zwall_45-zslit_45, rwall_45-rslit_45))
            
            
            rwall_195 = array([1.048,1.163,1.240,1.291,1.406,1.547,1.636,1.713,1.879,1.994,2.135,2.275 ]);
            zwall_195 = array([1.129,1.103,1.141,1.219,1.180,1.064,1.012,1.025,1.012,0.999,0.845,0.509 ]);
            zslit_195 = -0.834
            rslit_195 = 2.378
            
            xangle_195 = rad2deg(arctan2(zwall_195-zslit_195, rwall_195-rslit_195))
            
            
                       
            rxray  = r_[rslit_45, rslit_195*ones(12), rslit_195*ones(12)]*100
            zxray  = r_[zslit_45, zslit_195*ones(12), zslit_195*ones(12)]*100
            xangle = r_[xangle_45,xangle_195,xangle_195]

            
 
            #plot(vstack((rslit_195*ones(12),rwall_195)),vstack((zslit_195*ones(12),zwall_195)) )
            #Z = 0.000,0.964,0.968,1.001,1.019,1.077,1.070,1.096,1.113,1.138,1.147,1.165,1.217,\
            #1.217,1.162,1.162,1.163,1.165,1.166,1.166,1.169,1.172,1.176,1.183,1.183,1.185,\
            #1.188,1.191,1.196,1.202,1.208,1.214,1.221,1.231,1.238,1.244,1.254,1.278,1.290,\
            #1.331,1.347,1.348,1.348,1.348,1.348,1.348,1.348,1.310,1.310,1.292,1.095,1.077,\
            #1.077,1.040,0.993,0.709,0.519,0.389,0.400,0.222,0.133,0.044,-0.044,-0.133,-0.222,\
            #-0.400,-0.389,-0.973,-1.174,-1.211,-1.250,-1.250,-1.250,-1.329,-1.329,-1.363,\
            #-1.363,-1.363,-1.223,-1.223,-0.830,-0.800,-0.415,-0.400,-0.001,0.000

            #R= 1.016,1.016,1.016,1.016,1.016,1.016,1.016,1.016,1.016,1.016,1.016,1.012,1.001,\
            #1.029,1.042,1.046,1.056,1.097,1.108,1.116,1.134,1.148,1.162,1.181,1.182,1.185,1.190,\
            #1.195,1.201,1.209,1.215,1.222,1.228,1.234,1.239,1.242,1.248,1.258,1.263,1.280,1.280,\
            #1.280,1.310,1.328,1.361,1.380,1.419,1.419,1.372,1.372,1.608,1.647,1.785,2.070,2.128,\
            #2.245,2.323,2.377,2.362,2.364,2.365,2.365,2.365,2.365,2.364,2.362,2.377,2.134,1.786,\
            #1.768,1.768,1.682,1.372,1.372,1.420,1.420,1.273,1.153,1.016,1.016,1.016,1.016,1.016,\
            #1.016,1.016,1.016

            #figure(figsize=(4,7))
            #plot(c_[rxray,rxray2].T[:,:16], c_[zxray,zxray2].T[:,:16],'r')
            #plot(c_[rxray,rxray2].T[:,16:32], c_[zxray,zxray2].T[:,16:32],'y')
            #plot(c_[rxray,rxray2].T[:,32:48], c_[zxray,zxray2].T[:,32:48],'b')
            ##plot(c_[rxray,rxray2].T[:,48:64], c_[zxray,zxray2].T[:,48:64],'k')
            #ylim(-1.4,1.4)

            #plot(R,Z,'k')
            #axis('equal')
            #plt.axis([0.9, 2.45, -1.4, 1.4])

            #show()
                          #sxrgeom_plot3: load in geometry for 2017 SXR45U view chords
            zxray_new= array(( 73.1,)*8+( 78.1,)*4)
            rxray_new= array((230.4,)*8+(227.2,)*4)
            
            Rfin = array([2.19 , 1.98 ,  1.53, 1.12  , 1.01 , 1.01 , 1.01 ,1.01, 1.01, 1.01, 1.01, 1.01]);
            Zfin = array([-0.83 ,-1.06 ,-1.24 ,-1.34 ,-1.01 ,-0.65 ,-0.28 ,0.06, 0.20, 0.37, 0.61, 0.81]);

            xangle_new = rad2deg(arctan2(Zfin-zxray_new/100, Rfin-rxray_new/100))
       
            rxray_old = zeros(12) + 233.61
            zxray_old = zeros(12) + 78.62 
            xangle_old = array([ 227.94, 231.91, 236.02, 240.03, 244.11, 246.11, 
                248.14, 250.07, 252.04, 253.95, 255.83, 259.47 ])
            
            rxray  = r_[rxray_new, rxray_old, rxray_old]
            zxray  = r_[zxray_new, -zxray_old, -zxray_old]
            xangle = r_[xangle_new,-xangle_old,-xangle_old]
      
                    
            #print 'xxx', zxray
        elif shot >= 168847:
            
            #BUG it is wrong!!!!

            #sxrgeom_plot3: load in geometry for 2017 SXR45U view chords
            zxray_new= array(( 73.1,)*8+( 78.1,)*4)
            rxray_new= array((230.4,)*8+(227.2,)*4)
            
            Rfin = array([2.19 , 1.98 ,  1.53, 1.12  , 1.01 , 1.01 , 1.01 ,1.01, 1.01, 1.01, 1.01, 1.01]);
            Zfin = array([-0.83 ,-1.06 ,-1.24 ,-1.34 ,-1.01 ,-0.65 ,-0.28 ,0.06, 0.20, 0.37, 0.61, 0.81]);

            xangle_new = rad2deg(arctan2(Zfin-zxray_new/100, Rfin-rxray_new/100))
       
            rxray_old = zeros(12) + 233.61
            zxray_old = zeros(12) + 78.62 
            xangle_old = [ 227.94, 231.91, 236.02, 240.03, 244.11, 246.11, 
                248.14, 250.07, 252.04, 253.95, 255.83, 259.47 ]
            
            rxray  = r_[rxray_new, rxray_old, rxray_old]
            zxray  = r_[zxray_new, zxray_old, zxray_old]
            xangle = r_[xangle_new,xangle_old,xangle_old]

        elif shot >= 134917:
            rxray = zeros(12) + 233.61
            zxray = zeros(12) + 78.62 
            xangle = [ 227.94, 231.91, 236.02, 240.03, 244.11, 246.11, 
                248.14, 250.07, 252.04, 253.95, 255.83, 259.47 ]
        elif shot >= 130868:
            rxray = zeros(12) + 233.61
            zxray = zeros(12) + 78.62
            xangle=[236.124, 240.224, 241.936, 243.979, 246.014, 248.345, 250.409,
                252.682, 256.106, 259.752, 263.364, 267.078]
        elif shot >= 112250:
            rxray = zeros(12) + 233.61
            zxray = zeros(12) + 78.62
            xangle = [227.381,229.443,231.505,233.617,237.791,242.040,
                246.281,250.455,254.484,258.415,262.154,265.926]
        elif shot >= 108690:
            rxray = zeros(12) + 233.61
            zxray = zeros(12) + 78.62
            xangle = [227.401,229.523,231.579,233.633,237.713,241.750,
                245.905,249.886,253.810,257.618,261.196,264.833]
        elif shot >= 91000:
            rxray = zeros(12) + 232.0   
            zxray = zeros(12) + 78.56
        #     xangle = [234.60,236.60,238.60,240.70,244.70,248.70,
        #		252.60,256.40,260.10,263.60,266.40,270.00]
            xangle = [225.60,227.60,229.60,231.70,235.70,239.70,
                243.60,247.40,251.10,254.60,257.40,261.00]
        else:
            raise Exception('Toroidal X-ray is not available for shots earlier than 91000.')
        
        if len(xangle) == 12:
            xangle = r_[(xangle,)*3]
            rxray  = r_[(rxray, )*3]
            zxray  = r_[(zxray, )*3]

    else:
        #------------------------------------------------------------------------------
        # Poloidal SXR
        #------------------------------------------------------------------------------
        if shot < 80744:

            #  Old Xray arrays (0, -2 and 2)
            #  Print, 'Xraypaths: Old arrays for shots < 80744'
            xangle=[    120.,124.,128.,132.,136.,140.,144.,148.,  
                152.,156.,160.,164.,168.,172.,176.,180.,180., 
                184.,188.,192.,196.,200.,204.,208.,212.,216., 
                220.,224.,228.,232.,236.,240., 
                294.5,290.5,286.5,282.5,278.5,274.5,270.5,266.5, 
                262.5,258.5,254.5,250.5,246.5,242.5,238.5,234.5, 
                234.5,230.5,226.5,222.5,218.5,214.5,210.5,206.5, 
                202.5,198.5,194.5,190.5,186.5,182.5,178.5,174.5]

            zxray =[ -10.7,-10.7,-10.7,-10.7,-10.7,-10.7, 
                -10.7,-10.7,-10.7,-10.7,-10.7,-10.7,-10.7,-10.7, 
                -10.7,-10.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7, 
                -14.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7, 
                -14.7,-14.7, 
                130.1,130.1,130.1,130.1,130.1,130.1,130.1,130.1, 
                130.1,130.1,130.1,130.1,130.1,130.1,130.1,130.1, 
                132.6,132.6,132.6,132.6,132.6,132.6,132.6,132.6, 
                132.6,132.6,132.6,132.6,132.6,132.6,132.6,132.6]

            rxray =[ 248.9,248.9,248.9,248.9,248.9,248.9,248.9, 
                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, 
                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, 
                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, 
                197.6,197.6,197.6,197.6,197.6,197.6,197.6,197.6, 
                197.6,197.6,197.6,197.6,197.6,197.6,197.6,197.6,  
                194.5,194.5,194.5,194.5,194.5,194.5,194.5,194.5, 
                194.5,194.5,194.5,194.5,194.5,194.5,194.5,194.5] 

        elif shot < 91378:
            # sxr0s1-16, sxr0s17-32, sxr2s1-16, sxr2s17-32, DIII-D shot > 80744
            # 94/07/22 from Snider,modified for droop 1-6-95
            # took from EFIT /link/efit/sxr94.dat
            # print,'Xraypaths: for shots between 80744 and 91378.'


            r_0s1 = zeros(16) + 248.9
            z_0s1 = zeros(16) - 10.7
            theta_0s1 = [123.69,126.68,131.24,135.83,140.40,144.89,147.82,154.86,
                    157.87,160.92,163.98,167.06,170.12,173.15,176.13,179.05]

            r_0s17 = zeros(16) + 248.9
            z_0s17 = zeros(16) - 14.7
            theta_0s17= [180.33,183.22,186.18,189.19,192.24,195.30,198.36,201.40,
                    204.41,210.46,214.38,218.87,223.45,228.06,232.64,235.64]

            r_2s1 = zeros(16) + 197.6
            z_2s1 = zeros(16) + 130.1
            theta_2s1 = [290.17,285.65,281.05,277.98,274.92,271.89,268.91,265.98,
                    258.94,255.94,251.89,249.83,246.76,243.72,239.22,234.85]

            r_2s17 = zeros(16) + 194.5
            z_2s17 = zeros(16) + 132.6
            theta_2s17= [233.85,232.41,230.94,227.96,224.94,221.89,218.81,215.74,
                    212.70,209.68,202.71,199.71,196.72,193.69,190.63,187.56]

            xangle = r_[theta_0s1,theta_0s17,theta_2s1,theta_2s17]
            rxray  = r_[    r_0s1,    r_0s17,    r_2s1,    r_2s17]
            zxray  = r_[    z_0s1,    z_0s17,    z_2s1,    z_2s17]

        elif shot < 127721:
        #   New arrays (P1 and M1) from Robin Snider.   tbt

        #   Ted,  Below are the paths for the new x-ray arrays.  The data theta_side
        #   cooresponds to the m1 array.  i.e. the first element cooresponds to
        #   sx90rm1s1 and the last to sx90rm1s32.  The top array cooresponds to
        #   sx90rp1s1 ect.  The shot of interest is 89364.
        #
        #   The shot of interest should be 9137config8 instead of 89364.  The transition
        #   from old soft x-ray system and point names to the new ones for the 
        #   poloidal array was made on April,'97 after the RDP vent. (Snider 7-31-98)

        #   Print, 'Xraypaths: NEW arrays for shots >= 91378'

            theta_side = [98.8,101.6,104.6,107.5,110.5,113.5,115.      
                    ,116.5,118.,119.5,120.9,122.4,123.8,126.7,129.,131.9, 
                    134.8,139.2,143.7,148.1,152.5,156.8,158.8,163.,168.9 
                    ,174.8,180.7,185.,188.8,194.5,199.,207.9]

            z_side = [-77.65,-77.65,-77.65,-77.65,-77.65,-77.65,   
                    -77.65,-77.65,-77.65,-77.65,-77.65,-77.65,  
                    -77.65,-77.65,-77.65,-77.65,-77.65,-77.65,  
                    -77.65,-77.65,-77.65,-77.65,-81.36,-81.36,-81.36  
                    ,-81.36,-81.36,-81.36,-81.36,-81.36,-81.36,-81.36]

            r_side = [237.1,237.1,237.1,237.1,237.1,237.1,   
                    237.1,237.1,237.1,237.1,237.1,237.1,    
                    237.1,237.1,237.1,237.1,237.1,237.1,    
                    237.1,237.1,237.1,237.1,235.6,235.6,    
                    235.6,235.6,235.6,235.6,235.6,235.6     
                    ,235.6,235.6]

            theta_top = [260.9,258.,255.1,252.1,249.1,246.2,243.2   
                    ,240.2,237.3,234.5,230.7,227.9,223.5,219.1,214.6,  
                    210.2,205.9,203.1,201.4,198.6,195.7,192.8,189.9   
                    ,186.9,184.,179.6,175.3,171.5,168.7,165.8,162.9   
                    ,160.]

            z_top = [78.6,78.6,78.6,78.6,78.6,    
                    78.6,78.6,78.6,78.6,78.6,78.6,78.6,78.6,78.6,78.6  
                    ,78.6,78.6,78.6,82.3,82.3,82.3,82.3,82.3   
                    ,82.3,82.3,82.3,82.3,   
                    82.3,82.3,82.3,82.3,82.3]

            r_top = [236.9,236.9,236.9,236.9,236.9,236.9,   
                    236.9,236.9,236.9,236.9,236.9,236.9,  
                    236.9,236.9,236.9,236.9,236.9,236.9,  
                    235.3,235.3,235.3,235.3,235.3,235.3   
                    ,235.3,235.3,235.3,235.3,235.3,235.3  
                    ,235.3,235.3]

            xangle = r_[theta_top, theta_side]
            rxray  = r_[   r_top,     r_side]
            zxray  = r_[    z_top,     z_side]
            
            
     
        elif shot!=1e6 :
        # New poloidal SXR geometry for 2007 (E. Hollmann). xangle is defined
        # counterclockwise starting from R axis. Chords are RP1 (1:32), then RM1 (1:32).
        
            xangle = [270.8,267.2,263.3,259.2,254.8,250.3,245.6,240.9, 
                236.1,231.4,226.7,222.3,218.0,214.0,210.2,206.7, 
                217.8,214.3,210.6,206.5,202.3,197.8,193.2,188.4, 
                183.6,178.9,174.2,169.7,165.3,161.2,157.4,153.3, 
                89.2,92.8,96.7,100.8,105.2,109.7,114.4,119.1, 
                123.9,128.6,133.3,137.7,142.0,146.0,149.8,153.3, 
                142.2,145.7,149.4,153.5,157.7,162.2,166.8,171.6, 
                176.3,181.1,185.8,190.3,194.6,198.8,202.6,206.7]

            rxray = [230.4,230.4,230.4,230.4,230.4,230.4,230.4,230.4, 	
                230.4,230.4,230.4,230.4,230.4,230.4,230.4,230.4, 
                227.2,227.2,227.2,227.2,227.2,227.2,227.2,227.2, 	
                227.2,227.2,227.2,227.2,227.2,227.2,227.2,227.2, 
                230.4,230.4,230.4,230.4,230.4,230.4,230.4,230.4, 
                230.4,230.4,230.4,230.4,230.4,230.4,230.4,230.4, 
                227.2,227.2,227.2,227.2,227.2,227.2,227.2,227.2, 	
                227.2,227.2,227.2,227.2,227.2,227.2,227.2,227.2]
                
            zxray = [73.1,73.1,73.1,73.1,73.1,73.1,73.1,73.1, 
                73.1,73.1,73.1,73.1,73.1,73.1,73.1,73.1, 
                78.1,78.1,78.1,78.1,78.1,78.1,78.1,78.1, 
                78.1,78.1,78.1,78.1,78.1,78.1,78.1,78.1, 
                -73.1,-73.1,-73.1,-73.1,-73.1,-73.1,-73.1,-73.1, 
                -73.1,-73.1,-73.1,-73.1,-73.1,-73.1,-73.1,-73.1, 
                -78.1,-78.1,-78.1,-78.1,-78.1,-78.1,-78.1,-78.1, 
                -78.1,-78.1,-78.1,-78.1,-78.1,-78.1,-78.1,-78.1] 		


            xangle = array(xangle)
            #xangle[16:32],xangle[48:64] = copy(xangle[48:64][::-1]),copy(xangle[16:32][::-1])#NOTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! wrong geometry!!!
            
            
        else:
            
            
            xangle = [266.98,263.79,260.41,256.84,253.10,249.21,245.19,241.08,
                      236.92,232.75,228.63,224.58,220.64,216.86,213.24,209.80,
                      215.55,212.18,208.64,204.95,201.11,197.18,193.16,189.11,
                      185.07,181.07,177.15,173.35,169.69,166.19,162.86,159.72,
                      91.03,94.45,98.08,101.90,105.89,110.03,114.28,118.59,
                      122.92,127.21,131.43,135.52,139.46,143.21,146.77,150.12,
                      146.00,149.36,152.91,156.64,160.52,164.54,168.66,172.84,
                      177.02,181.18,185.27,189.24,193.07,196.74,200.22,203.50]
            
            
            rxray = array([230.4]*16+[227.2]*16+[230.4]*16+[227.2]*16)
            zxray = array([ 73.1]*16+[ 78.1]*16+[-73.1]*16+[-78.1]*16)

            ##; Bottom channels
            #rwallbot = array([ 221, 209, 197, 182, 163, 147, 129, 112, 
                        #102, 102, 102, 102, 102, 102, 102,102       ]) 
            #zwallbot = array([ -78, -99,-107,-115,-136,-136,-136,-133, 
                        #-118, -94, -73, -56, -40, -27, -14,-1       ]) 
            ##; Top channels
            #rverttop = array((227.2,)*16)
            #zverttop =  array((78.1,)*16)
            #rwalltop = array([ 102, 102, 102, 102, 102, 102, 102, 102, 102, 
                        #102, 102, 102, 102, 102, 117, 122       ]) 
            #zwalltop = array([ -2,  7,  16,  25,  35,  44,  53,  61,  71, 
                        #80,  88,  97, 106, 115, 119, 123       ]) 
            

            #rxray = r_[rvertbot, rverttop,rvertbot, rverttop  ]
            #zxray = r_[zvertbot, zverttop, -zvertbot, -zverttop ]

            #rxray2 = r_[rwallbot, rwalltop,rwallbot, rwalltop  ]
            #zxray2 = r_[zwallbot, zwalltop, -zwallbot, -zwalltop ]

            #xangle = rad2deg(arctan2(zxray2-zxray,rxray2-rxray))
    


    rxray = array(rxray)/100.
    zxray = array(zxray)/100.

    angle =  deg2rad(xangle) 
    
    
    rxray2 = rxray  + 3 * cos(angle)
    zxray2 = zxray  + 3 * sin(angle)
    
    #Z = 0.000,0.964,0.968,1.001,1.019,1.077,1.070,1.096,1.113,1.138,1.147,1.165,1.217,\
    #1.217,1.162,1.162,1.163,1.165,1.166,1.166,1.169,1.172,1.176,1.183,1.183,1.185,\
    #1.188,1.191,1.196,1.202,1.208,1.214,1.221,1.231,1.238,1.244,1.254,1.278,1.290,\
    #1.331,1.347,1.348,1.348,1.348,1.348,1.348,1.348,1.310,1.310,1.292,1.095,1.077,\
    #1.077,1.040,0.993,0.709,0.519,0.389,0.400,0.222,0.133,0.044,-0.044,-0.133,-0.222,\
    #-0.400,-0.389,-0.973,-1.174,-1.211,-1.250,-1.250,-1.250,-1.329,-1.329,-1.363,\
    #-1.363,-1.363,-1.223,-1.223,-0.830,-0.800,-0.415,-0.400,-0.001,0.000

    #R= 1.016,1.016,1.016,1.016,1.016,1.016,1.016,1.016,1.016,1.016,1.016,1.012,1.001,\
    #1.029,1.042,1.046,1.056,1.097,1.108,1.116,1.134,1.148,1.162,1.181,1.182,1.185,1.190,\
    #1.195,1.201,1.209,1.215,1.222,1.228,1.234,1.239,1.242,1.248,1.258,1.263,1.280,1.280,\
    #1.280,1.310,1.328,1.361,1.380,1.419,1.419,1.372,1.372,1.608,1.647,1.785,2.070,2.128,\
    #2.245,2.323,2.377,2.362,2.364,2.365,2.365,2.365,2.365,2.364,2.362,2.377,2.134,1.786,\
    #1.768,1.768,1.682,1.372,1.372,1.420,1.420,1.273,1.153,1.016,1.016,1.016,1.016,1.016,\
    #1.016,1.016,1.016

    #figure(figsize=(4,7))
    #plot(c_[rxray,rxray2].T[:,:16], c_[zxray,zxray2].T[:,:16],'r')
    #plot(c_[rxray,rxray2].T[:,16:32], c_[zxray,zxray2].T[:,16:32],'y')
    #plot(c_[rxray,rxray2].T[:,32:48], c_[zxray,zxray2].T[:,32:48],'b')
    #plot(c_[rxray,rxray2].T[:,48:64], c_[zxray,zxray2].T[:,48:64],'k')
    #ylim(-1.4,1.4)

    #plot(R,Z,'k')
    #axis('equal')
    #plt.axis([0.9, 2.45, -1.4, 1.4])

    #show()

    return rxray,zxray,rxray2,zxray2,angle




### FROM OMFIT !!!!!!
#-----------------------------
# Get filter settings by year
# From /u/hollmann/SXR/SXRsettingsPA.dat
# and /u/hollmann/matlab/sxr_calib.m
#-----------------------------

shot_year = (
    (0, 2007),
    (127330, 2007),
    (131060, 2008),
    (135130, 2009),
    (140840, 2010),
    (143710, 2011),
    (148160, 2012),
    (152130, 2013),
    (156200, 2014),
    (160940, 2015),
    (164780, 2016),
    (168440, 2017),
    (174580, 2018))

def get_filter(shot, PA):
    
    #find a year corresponding to the discharge number, BUG use SQL? 
    year = 0
    for shot_, year_ in shot_year:
        if  shot < shot_: break
        year = year_
        
    #print year
    
    #exit()
        
    if PA and year < 2007:
        raise Exception('No SXR filter information!')
        #OMFITx.End()
        
    if PA and year in (2007,2008):
        settings=[{'description':'Closed', 'aperture':'closed', 'foil':'none', 'filter':'none','pinhole diameter':0.0},
                  {'description':'Disruption prad', 'aperture':'200 um pinhole', 'foil':'50 um thick moly', 'filter':'none','pinhole diameter':200.0},
                  {'description':'ELM prad', 'aperture':'300 um x 2 mm slit', 'foil':'50 um thick moly', 'filter':'none','pinhole diameter':870.0},
                  {'description':'Low Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'50 um thick moly', 'filter':'12 um diamond','pinhole diameter':1.95e3},
                  {'description':'Mid Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'50 um thick moly', 'filter':'15 um diamond','pinhole diameter':1.95e3},
                  {'description':'Hi Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'50 um thick moly', 'filter':'32 um diamond','pinhole diameter':1.95e3}]
    if PA and year in (2009,2010):
        settings=[{'description':'Closed', 'aperture':'closed', 'foil':'none', 'filter':'none','pinhole diameter':0.0},
                  {'description':'Disruption prad', 'aperture':'200 um pinhole', 'foil':'?', 'filter':'none','pinhole diameter':200.0},
                  {'description':'ELM prad', 'aperture':'300 um x 2 mm slit', 'foil':'?', 'filter':'none','pinhole diameter':870.0},
                  {'description':'Low Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'?', 'filter':'12 um diamond','pinhole diameter':1.95e3},
                  {'description':'Mid Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'?', 'filter':'18 um diamond','pinhole diameter':1.95e3},
                  {'description':'Hi Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'?', 'filter':'10 um Ti','pinhole diameter':1.95e3}]
    if PA and year == 2011:
        settings=[{'description':'Closed', 'aperture':'closed', 'foil':'none', 'filter':'none','pinhole diameter':0.0},
                  {'description':'Disruption prad', 'aperture':'200 um pinhole', 'foil':'?', 'filter':'none','pinhole diameter':200.0},
                  {'description':'ELM prad', 'aperture':'300 um x 2 mm slit', 'foil':'?', 'filter':'none','pinhole diameter':870.0},
                  {'description':'Low Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'?', 'filter':'12 um Be','pinhole diameter':1.95e3},
                  {'description':'Mid Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'?', 'filter':'25 um Be','pinhole diameter':1.95e3},
                  {'description':'Hi Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'?', 'filter':'125 um Be','pinhole diameter':1.95e3}]
    if PA and year in (2012,2013):
        settings=[{'description':'Closed', 'aperture':'closed', 'foil':'none', 'filter':'none','pinhole diameter':0.1},
                  {'description':'Disruption prad', 'aperture':'200 um pinhole', 'foil':'50 um thick SS', 'filter':'none','pinhole diameter':200.0},
                  {'description':'ELM prad', 'aperture':'300 um x 2 mm slit', 'foil':'25 um thick moly', 'filter':'none','pinhole diameter':870.0},
                  {'description':'Low Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'1.3 mm thick tantalum', 'filter':'12 um Be','pinhole diameter':1.95e3},
                  {'description':'Visible light', 'aperture':'1 mm x 3 mm slit', 'foil':'25 um thick moly', 'filter':'1 mm fused silica','pinhole diameter':1.95e3},
                  {'description':'Hi Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'25 um thick moly', 'filter':'125 um Be','pinhole diameter':1.95e3}]
    if PA and year in r_[2014:2019]:
        settings=[{'description':'Closed', 'aperture':'closed', 'foil':'none', 'filter':'none','pinhole diameter':0.0},
                  {'description':'Disruption prad', 'aperture':'200 um pinhole', 'foil':'50 um thick SS', 'filter':'none','pinhole diameter':200.0},
                  {'description':'ELM prad', 'aperture':'300 um x 3 mm slit', 'foil':'25 um thick moly', 'filter':'none','pinhole diameter':1070.},#870.0},
                  {'description':'Low Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'1.3 mm thick tantalum', 'filter':'12 um Be','pinhole diameter':1.95e3},
                  {'description':'Intermediate prad', 'aperture':'400 um pinhole', 'foil':'50 um thick SS', 'filter':'none','pinhole diameter':400.0},
                  {'description':'Hi Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'25 um thick moly', 'filter':'125 um Be','pinhole diameter':1.95e3}]
        
        
        
    if not PA and shot < 168847:
        settings=[{'description':'140um Be', 'aperture':' '     ,'foil':'none',     'filter':'140um Be','pinhole diameter':3.91e3}]

    elif not PA and shot < 172691:
        settings=[
                {'description':'Closed',   'aperture':'closed','foil':'none',     'filter':'none','pinhole diameter':0.0},
                {'description':'10 um Ti', 'aperture':'1 mm x 1.2 mm slit '     ,'foil':'10 um Ti', 'filter':'none','pinhole diameter':3.91e3},
                {'description':'ND1',      'aperture':'1 mm x 1.2 mm slit '     ,'foil':'none',     'filter':'none','pinhole diameter':3.91e3},
                {'description':'ND2',      'aperture':'1 mm x 1.2 mm slit '     ,'foil':'none',     'filter':'none','pinhole diameter':3.91e3},
                {'description':'25um Be',  'aperture':'1 mm x 1.2 mm slit '     ,'foil':'none',     'filter':'25um Be','pinhole diameter':3.91e3},
                {'description':'125um Be', 'aperture':'1 mm x 1.2 mm slit '     ,'foil':'none',     'filter':'125um Be','pinhole diameter':3.91e3}]

        
    elif not PA:
        settings=[
                {'description':'Closed',          'aperture':'closed','foil':'none',    'filter':'none',     'pinhole diameter':0.0},
                {'description':'dummy for bakes', 'aperture':'1 mm x 1.2 mm slit ',     'foil':'none',    'filter':'100 um SS','pinhole diameter':3.91e3},
                {'description':'default SXR',     'aperture':'1 mm x 1.2 mm slit ',     'foil':'none',    'filter':'125 um Be','pinhole diameter':3.91e3},
                {'description':'low Te SXR',      'aperture':'1 mm x 1.2 mm slit ',     'foil':'none',    'filter':'25um Be',  'pinhole diameter':3.91e3},
                {'description':'disruptions Prad','aperture':'1 mm x 1.2 mm slit ',     'foil':'none',    'filter':'none',     'pinhole diameter':3.91e3},
                {'description':'ELM Prad',        'aperture':'1 mm x 1.2 mm slit ',     'foil':'none',    'filter':'none',     'pinhole diameter':3.91e3}]

                

    return settings

#-----------------------------
# SXR Calibration to turn Volts into W/cm**2
#-----------------------------
def get_calib(shot,calib_path,cam):
    # Get the filter settings options for our shot
    if verbose: print('Getting SXR Filter Settings')
    
    PA = cam in ('90RM1','90RP1')
    filter_settings = get_filter(shot,PA)
    # Get electronic settings and filter setting

    if PA:
        SXRsettings = calib_path+os.sep+'SXRsettingsPA.dat'  
    else:
        SXRsettings = calib_path+os.sep+'SXRsettings45U.dat'  #BUG!!!
    index = open(SXRsettings,'r')

    shots, Rc, Gain,Filt = [],[],[],[]
    do_read = False
    for line in index:
        if do_read:
            P = cam == '90RP1'
            line = line.split()
            shots.append(int(line[0]))
            
            Rc.append(float(line[1+P].replace('k','e3')))
            Gain.append(float(line[2+PA+P]))
            Filt.append(int(line[3+PA*2+P]))
        if line[:4] == 'shot':
            do_read = True
    index.close()
    # Find settings for our shot
    wh = where(array(shots) >= shot)[0]
    # Go back one for the change shot
    try:
        whs = wh[0]-1
    except Exception:
        print('Calibration shot is last shot!')
        whs = len(shots)-1

    Rc, Gain, filter = Rc[whs], Gain[whs], Filt[whs]
    

    if (shot < 168847 and cam == '45R1') or cam in ('165R1','195R1'):
        #old TA diagnostics
        Rc, Gain, pinhole,filt   = 500e3, 1,  3.91e3,1
    else:
        # Get pinhole diameter (um)
        pinhole = filter_settings[filter]['pinhole diameter']


    if pinhole == 0:
        print('Pinhole closed '+cam)
        pinhole = 1.95e3

    print('%s\tRc:%.0e  Gain:%d pinhole:%.2g filter:%d'%(cam, Rc, Gain, pinhole,filter))


    return Rc, Gain, pinhole, filter

def get_calib_fact(shot, geometry_path,  toroidal=False):
    ## --------------------------
    # Get subsystem calibrations
    # --------------------------
    if verbose: print('Getting SXR Calibration Settings')
    calib_dict = {}

    if toroidal:
        cams = '45R1','165R1','195R1'
    else:
        cams = '90RM1','90RP1'


    for cam in cams:
        #BUG from for Tor cameras!!!!!!
        calib_path = '/fusion/projects/diagnostics/sxr/'

        try:
            resistor, gain, pinhole, filter = get_calib(shot,calib_path,cam)
        except:
            #if cam == '195R1': shot = min(shot, 160000) #to prevent loading of calibration for new U45 camera
            resistor, gain, pinhole, filter = get_calib(shot,geometry_path,cam)

        # responsivity of AXUV photodiodes [A/W] from E Hollman
        eff_resp =  0.27
        
        # 50 ohm termination divides by 2
        term = 0.5  
        
        if cam in ('90RM1','90RP1'): #poloidal array
            pol_array = True
            nch = 16
            # Effective pinhole thickness (microns)
            pinhole_thick = 25.#50.0
            # Distance center to pinhole (cm)
            ctr_to_pinhole = 2.95#3.00
            # Element area (m^2)
            el_area = 2 * 5/1e6
            # Aperture area (m^2)
            ap_area = 0.25*pi*(pinhole/1e6)**2
            # Distance from center (cm) in eight steps
            dist_from_ctr = ( arange(-nch/2,nch/2) + 0.5 ) * 0.212 #   % [cm]
            tanpsi = dist_from_ctr / ctr_to_pinhole
            cos4 = (tanpsi**2 + 1.)**(-2.)
            thick_factor = abs(tanpsi) * (-4./pi) * (pinhole_thick/pinhole) + 1.

            #measured for filter 5 (high Te)
            if cam == '90RM1':
                etendue  = array( [2.141,2.338,2.528,2.702,2.882,2.695,2.728,2.985,
                                   2.989,2.971,2.858,2.711,2.538,2.349,2.152,1.900,
                                   2.074,2.331,2.589,2.836,3.056,3.234,3.356,3.41,
                                   3.382,3.282,3.122,2.914,2.675,2.419,2.161,1.91])*1e-8
            if cam == '90RP1':
                etendue  = array( [1.845,2.073,2.302,2.523,2.72,2.671,2.813,3.015,
                                   3.031,2.948,2.811,2.631,2.421,2.195,1.965,1.74,
                                   1.881,2.182,2.425,2.663,2.882,3.069,3.211,3.296,
                                   3.314,3.259,3.142,2.975,2.769,2.539,2.297,2.055])*1e-8
                
            if filter == 3:
                #measured, thick slits for low Te
                if cam == '90RM1':
                    etendue  = array( [1.749,2.18,2.442,2.646,2.754,2.558,2.592,2.783,
                                       2.824,2.833,2.783,2.646,2.504,2.309,2.031,1.765,
                                       2.347,2.646,2.895,3.098,3.244,3.352,3.414,3.447,
                                       3.447,3.352,3.211,3.053,2.866,2.625,2.23,1.661])*1e-8
                if cam == '90RP1':
                    etendue  = array( [1.628,2.068,2.356,2.546,2.648,2.607,2.713,2.847,
                                       2.945,2.859,2.758,2.632,2.498,2.315,2.067,1.645,
                                       1.877,2.076,2.292,2.479,2.658,2.766,2.865,2.931,
                                       2.946,2.905,2.803,2.654,2.492,2.23,1.873,1.279])*1e-8
       
            else:
                #correction for different area of the slit 
                etendue *= pinhole**2/1.95e3**2
       
                
            
        if cam in ('45R1','165R1','195R1'): #toroidal array
            pol_array = False
            nch = 20
            # Effective pinhole thickness (microns)
            pinhole_thick = 0.
            # Distance center to pinhole (cm)
            ctr_to_pinhole = 2.7
            # Element area (m^2)
            el_area = 4.1*0.75/1e6
            # Aperture area (m^2)
            ap_area = 0.25*pi*(pinhole/1e6)**2# 12*1/1e6
            #pinhole = ap_area*
            # Distance from center (cm) in eight steps
            dist_from_ctr = ( arange(-nch//2,nch//2) + 0.5 ) * 0.095 #   % [cm]
            tanpsi = dist_from_ctr / ctr_to_pinhole
            cos4 = (tanpsi**2 + 1.)**(-2.)
            thick_factor = abs(tanpsi) * (-4./pi) * (pinhole_thick/pinhole) + 1.

            # overall factor to account for optics...units: m-2 sr-1  it is 1/etendue
            etendue = cos4*thick_factor*el_area*ap_area/ctr_to_pinhole**2*1e4
          
        # Compute calibration
               
        # From V -> W/m^2
        
        #0.35 A/W    U = IR 
        #print 0.35*5e-4*(resistor*gain)
        

        calib = 4*pi/term/(resistor*gain)/eff_resp/etendue
        print('calib '+cam+' %.2e  W/m^2/V  '%mean(calib)+ ' %.2e  W/V  '%mean(calib*el_area))
        
        
        #W = C*plocha*Volty
        #volty = W/plocha/C
        
        #plot(calib);show()
        #calib[:] = 1
        #### Done computing calibration ####
        if pol_array:
            calib_dict[cam+'a'] = calib[:16]
            calib_dict[cam+'b'] = calib[16:]
        else:
            active = [18,17,15,13,11,10,9,8,7,6,5,3]  #from Eric Hollman

            calib_dict[cam] = calib[active]

        
    return  calib_dict

        

### UP TO HERE FROM OMFIT !!!!!!

def mds_load(params):
    #print TDI
    (mds_server, tree,shot,  TDI) = params
    MDSconn = mds.Connection(mds_server )
    output = []
    
    if tree is not None: 
        MDSconn.openTree(tree,shot)
        
    
    for tdi in TDI:
        try:
            output.append(MDSconn.get(tdi).data())
        except:
            print('missing data:'+ tdi)
            output.append(array(0, ndmin=1))

            
    if tree is not None: 
        MDSconn.closeTree(tree,shot)

    return output

from multiprocessing.pool import  Pool

def mds_par_load(mds_server,  TDI,  numTasks,tree=None,shot=None):

    #load a junks of a single vector

    TDI = array_split(TDI, min(numTasks, len(TDI)))

    args = [(mds_server, tree,shot, tdi) for tdi in TDI]

    pool =  Pool(len(args))
    

    out = pool.map(mds_load,args)
    pool.close()
    pool.join()

    output = []
    for o in out:output+= o
    
    return  output
    













class loader_SXR():


    def __init__(self,shot, geometry_path,MDSconn, fast_data, toroidal):
        
        self.shot = shot
        self.geometry_path = geometry_path
        self.MDSconn = MDSconn
        
        if self.MDSconn is None and fast_data:
            raise Exception('MDS connection is not available')
        
        #exit()
        self.fast_data = fast_data
        self.tree = 'spectroscopy'
        self.toroidal = toroidal


        self.filter_data = True
        self.calib = get_calib_fact(self.shot,geometry_path,toroidal )
        self.wrong_dets_damaged = []
        # Geometry from
        # /usc-data/c/idl/source/efitviewdiagnoses/DIII-D/xraypaths.pro
        # Start with R+1 #1-32 then R-1 #1-32
        # angle of view from R in degrees
      
        if not os.path.exists(self.geometry_path+'data'):
            os.mkdir(self.geometry_path+'data')
 
        r_p,z_p,r2_p,z2_p,xangle_p = xraypaths(1e6, toroidal=False)
        r_t,z_t,r2_t,z2_t,xangle_t = xraypaths(self.shot, toroidal=True)


        self.Phi = {'165R1':[165]*12, '195R1':[195]*12,
                    '45R1':[45]*12, '90RM1a':[90]*16,'90RM1b':[90]*16, '90RP1a':[90]*16,'90RP1b':[90]*16  }           
        self.R_start= {"90RM1a":r_p[32:48],"90RM1b":r_p[48:],
                       "90RP1a":r_p[:16],"90RP1b":r_p[16:32],
                       '45R1':r_t[0:12],'165R1':r_t[12:24],'195R1':r_t[24:36]}
        self.R_end=   {"90RM1a":r2_p[32:48],"90RM1b":r2_p[48:],
                       "90RP1a":r2_p[:16],"90RP1b":r2_p[16:32],
                       '45R1':r2_t[0:12],'165R1':r2_t[12:24],'195R1':r2_t[24:36]}
        self.z_start= {"90RM1a":z_p[32:48],"90RM1b":z_p[48:],
                       "90RP1a":z_p[:16],"90RP1b":z_p[16:32],
                       '45R1':z_t[0:12],'165R1':z_t[12:24],'195R1':z_t[24:36]}
        self.z_end=   {"90RM1a":z2_p[32:48],"90RM1b":z2_p[48:],
                       "90RP1a":z2_p[:16],"90RP1b":z2_p[16:32],
                       '45R1':z2_t[0:12],'165R1':z2_t[12:24],'195R1':z2_t[24:36]}
        self.theta=   {"90RM1a":xangle_p[32:48],"90RM1b":xangle_p[48:],
                       "90RP1a":xangle_p[:16],"90RP1b":xangle_p[16:32],
                       '45R1':xangle_t[0:12],'165R1':xangle_t[12:24],'195R1':xangle_t[24:36]}
        
        self.cam_geom = {'pol':{'Dx': 2.,
                                'Dy': 5.,
                                'Px': 1.,
                                'Py': 3.,
                                'Delta': 29.},
                        'tor':{ 'Dx': 0.75, 
                                'Dy': 4.1,
                                'Px': 1., 
                                'Py': 12, 
                                'Delta': 27.}}

        self.wrong_dets = []
        self.detectors_dict=OrderedDict()

        if toroidal:
            cams = '45R1','195R1',
            self.n_diods = 12
            n_chunks = 5
            for c in cams: 
                self.detectors_dict[c] = []
                for j in range(self.n_diods):
                    self.detectors_dict[c].append(c+'_'+str(j+1))
            self.calb_0 = 1,1  
   
                    
        else:
            cams = '90RP1a','90RP1b','90RM1a','90RM1b'

            self.n_diods = 16
            n_chunks = 4
            for i,c in enumerate(cams):
                self.detectors_dict[c] = []
                for j in range(self.n_diods):
                    self.detectors_dict[c].append(c+'_'+str(self.n_diods*(i%2)+j+1))
            self.calb_0 = 1,1,1,1 

        self.cam_ind = OrderedDict()
        count = 0
        for i, c in enumerate(cams):
            self.cam_ind[c] = list(range(count,count+self.n_diods))
            count+= self.n_diods
        
        self.all_los = hstack([self.detectors_dict[c] for c in cams])
        self.n_dets = len(cams)
        self.nl = self.n_diods*self.n_dets
        self.dets_index = [array(v) for k,v in self.cam_ind.items()]
        self.dets = arange(self.nl)
        
        self.fast_downsample = 1
        self.fast_data_downsapled=False
   
        if fast_data:
            MDSserver = self.MDSconn.hostspec

            suffix = 'TA' if toroidal else 'PA'
            catch_path = os.path.join(self.geometry_path,'data','%d_%s_fast.npz'%(self.shot, suffix))
            if config.useCache:
                try:
                    data_file = load(catch_path,allow_pickle=True,encoding='latin1' )
                    self.tvec_fast = data_file['tvec'] 
                except:
                    pass
        
            if not hasattr(self, 'tvec_fast'):
                print('fetching time vec')
                #exit()

                cam = cams[0][:-2] if toroidal else '90'+cams[0][3]
                ch = 1
                
                if self.shot < 136162:   # old fast data
                    TDI = ('dim_of(PTDATA2("SX%sF%.2d",%d))'%(cam,ch,self.shot),)
                    self.tvec_fast = mds_par_load(MDSserver, TDI,1)
                else:# new fast data (DTACQ system)
                    
                    TDI = ['dim_of(PTDATA2("SX%sF%.2d_%d",%d))'%(cam,ch,i,self.shot) for i in range(n_chunks)]
                    self.tvec_fast = mds_par_load(MDSserver, TDI,8)
                    if len(self.tvec_fast[0]) == 1:
                        #print TDI
                        raise Exception("fast data are not availible!")
    
                
                for i in range(n_chunks): 
                    self.tvec_fast[i] /= 1e3 #s
                    if self.fast_downsample > 1:
                        self.tvec_fast[i] = self.tvec_fast[i].reshape(-1,self.fast_downsample).mean(1) 

                
            
            self.tvec = hstack(self.tvec_fast)#.reshape(-1,self.fast_downsample).mean(1) #dowsample
        else:
            suffix = 'TA' if self.toroidal else 'PA'
            catch_path = os.path.join(self.geometry_path,'data','%d_%s.npz'%(self.shot, suffix))
            if os.path.isfile(catch_path) and config.useCache:
                self.tvec = load(catch_path,allow_pickle=True,encoding='latin1')['tvec']
                
            elif self.MDSconn is not None:
                cam = cams[0] if toroidal else cams[0][:-1]
                try:
                    #print('DO not load fast dowsampled data')
                    #eee
                    #downsampled fast data, if availible
                    self.MDSconn.openTree('spectroscopy', self.shot)
                    self.tvec = MDSconn.get('dim_of(\\SX165R1F:SX165R1F05)').data()/1e3
                    self.fast_data_downsapled=True
                    print('===============  fast downsampled data!! =======')

                except:
                    ##slow diag 
                
                    TDIcall = 'PTDATA2("SX%sS%s",%d)'%(cam,1,self.shot)
                    self.tvec = self.MDSconn.get('dim_of('+TDIcall+')').data()/1e3

                    if len(self.tvec) == 1:
                        self.cache= None
                        raise Exception("slow data are not availible!")
            else:
                raise Exception("MDS is not connected")

 

    def load_geom(self,path):
        
        #separate equlibrium for each campaign 
        ##load corrections in degrees of the camera position
        self.geometry_version = 1 if self.shot > 168847 else 0 
        
        suffix = 'tor' if self.toroidal else 'pol'
        corrections = 'det_pos_corr_null'
        corrections = 'det_pos_corr_new'

        pos_corr = loadtxt(path+'/'+corrections+'.txt',
                           dtype={'names': ('det', 'angle'),'formats': ('U6',  'f4')})
        pos_corr =  {k:item for k,item in pos_corr}
 
    
        coord_dict = OrderedDict()
 
        ###prepare files with cordinates
        
        for icam,(det,channels) in enumerate(self.detectors_dict.items()):
   
            xfile = open(self.geometry_path+'detector_%s_x.txt'%det,'w')
            yfile = open(self.geometry_path+'detector_%s_y.txt'%det,'w')
            dist_file = open(self.geometry_path+'dist_%s.txt'%det,'w')

            coord_dict[det] = []
            verts = []
  
            
            R1 = self.R_start[det]
            R2 = self.R_end[det]
            Z1 = self.z_start[det]
            Z2 = self.z_end[det]
            THETA = self.theta[det]

            m_alpha = mean(abs(abs(arctan2(Z2-Z1,R2-R1))-pi/2))

            for ich,ch in enumerate(channels):
                r1 = R1[ich]
                z1 = Z1[ich]
                r2 = R2[ich]
                z2 = Z2[ich]
                theta = THETA[ich]

                coord_dict[det].append([[r1,r2, z1,z2]])

                dAngle = deg2rad(theta-mean(THETA))
                if dAngle > pi: dAngle-= 2*pi
     

                alpha = arctan2(z2-z1,r2-r1)
                #print det
                if det  in pos_corr:
                    alpha+= deg2rad(pos_corr[det])
                #print(det, pos_corr)
                L = hypot(r2-r1, z2-z1)
                if L == 0: L = 1;alpha = 1
                
                cam_geom = self.cam_geom[suffix]
                #camera geometry P - pinhole, D - diode
                Delta = cam_geom['Delta']#self.geometry['Foc_Len'][sig]#14.#mm
                Px = cam_geom['Px']#self.geometry['P_Width'][sig] # 0.3#mm
                Py = cam_geom['Py']#self.geometry['D_Length'][sig] # 4.6#mm
                Dx = cam_geom['Dx']#self.geometry['D_Width'][sig] #  0.96 mm
                Dy = cam_geom['Dy']#self.geometry['P_Length'][sig] # 5.0#mm
                
      
                
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
            dist_file.close()
 

    
    def get_data_fast(self,tmin=-infty,tmax=infty,calib=True):
    
        num_MDS_Tasks = 8
        MDSserver = self.MDSconn.hostspec


        if self.shot < 136162:
            raise Exception('loading from old DAS is not implemented')
        indmin = where([t[-1] > tmin for t in self.tvec_fast])[0][0]
        indmax = where([t[ 0] < tmax for t in self.tvec_fast])[0][-1]+1
        index = unique(r_[0,indmin:indmax])
        
        #print tmin, tmax,index 
        
        cam_ind = {c:slice(ii[0],ii[-1]+1) for c,ii in self.cam_ind.items()} 
        
        tvec = hstack(self.tvec_fast[indmin:indmax])  
        
        #self.cache_fast
        
        suffix = 'TA' if self.toroidal else 'PA'
        catch_path = os.path.join(self.geometry_path,'data','%d_%s_fast.npz'%(self.shot, suffix))
        if not hasattr(self,'cache_fast') and config.useCache:
            if os.path.isfile(catch_path):
                
                try:
                    data_file = load(catch_path,allow_pickle=True,encoding='latin1')
                    self.cache_fast = data_file['cache_fast'].item()
                except Exception as e:
                    print(e)
        
            
        if not hasattr(self,'cache_fast') :
            self.cache_fast = {}
            
        #load only signals which are not cached. 
        TDI = []
        for i in index:
            if i not in self.cache_fast:
                for cam, channels in self.detectors_dict.items():
                    cam = cam[:-2] if self.toroidal else '90'+cam[3]
                    for ch in channels:
                        ch = int(ch.split('_')[-1])
                        TDIcall = 'PTDATA2("SX%sF%.2d_%d",%d)'%(cam,ch,i,self.shot)
                        TDI.append(TDIcall)
                        
        if len(TDI) > 0:
            print('fetching data'+ str(TDI))
            t = time.time()
            data = mds_par_load(MDSserver, TDI,  num_MDS_Tasks)
            print('done in %.2fs'%( time.time()-t))

        j = 0
        for i in index:
            if i not in self.cache_fast:
                self.cache_fast[i] = {c:{} for c in list(self.detectors_dict.keys())}
                for cam, channels in self.detectors_dict.items():
                    for ch in channels:#downsample to 125kHz
                        #supress spikes?
                        if self.fast_downsample  > 1:
                            self.cache_fast[i][cam][ch] = median(data[j].reshape(-1,self.fast_downsample),1)#.mean(1)
                        else:
                            self.cache_fast[i][cam][ch] = data[j]
                        #print self.cache_fast[i][cam][ch].dtype
                        j+= 1
        
        #print catch_path
        if len(TDI) > 0 and config.useCache:
            savez_compressed(catch_path, tvec=self.tvec_fast, cache_fast=self.cache_fast)
            #save the cache
            


        #merge data to a single array

        data = []
        for cam, channels in self.detectors_dict.items():
            for ch in channels:
               data.append(hstack([self.cache_fast[i][cam][ch] for i in range(indmin,indmax)]))
        
        data  = vstack(data).T
        data_err = zeros_like(data)
    

        #remove corrupted channels
        self.hardcoded_corrections(tvec, data,data_err)
    
        
        #estimate random noise
        try:
            nt = len(tvec)
            ind_t = slice(0,nt)
            if nt > 2000:  ind_t = r_[0, unique(randint(nt,size=1000)),nt-1]
            from shared_modules import fast_svd
            U,S,V = fast_svd(data[ind_t],min(30,nt//3))

            #assume that the differnce between data and SVD retrofit is only noise 
            data_err+= std(data[ind_t]-dot(U, V*S[:,None]),0)
        except:pass


        imin,imax = tvec.searchsorted([tmin,tmax])
        ind = slice(imin,imax)
        #if tmin > 1.5:
        ##print data.shape
            #d0 = data[(tvec > 1.8 )&(tvec < 1.81)].mean(0)[None]
            #d1 = data[(tvec > 1.99)&(tvec < 2.)].mean(0)[None]
            ##print d0
            ##print d1
            #data += (d1+d0)/2-(tvec[:,None]-1.8)/(2.-1.8)*(d1-d0)+d0

        #print data
        tvec, data, data_err = tvec[ind], data[ind], data_err[ind]

        #calibrate raw signals
        if calib:
            tvec_offset, data_offset,_ = self.get_data_fast(tmin=-infty,tmax=.2, calib=False)
            data -= data_offset.mean(0)
            data_offset-= data_offset.mean(0)
            data_err+= data_offset.std(0)/10  
            
            for cam, ind in cam_ind.items():
                if cam in self.calib:
                    data[:,ind]  *= self.calib[cam]
                    data_err[:,ind] *= self.calib[cam]
                else:
                    print('Warning: calibration for camera %s is not availible'%cam)

        return tvec, data, data_err

                 

 
    def get_data(self,tmin,tmax,fetch_errors=True):
        
        if self.fast_data:
            return self.get_data_fast(tmin,tmax)
     
        #store data for faster loading 
        suffix = 'TA' if self.toroidal else 'PA'
        catch_path = os.path.join(self.geometry_path,'data','%d_%s.npz'%(self.shot, suffix))

        if not hasattr(self,'cache') and config.useCache:
        
            if os.path.isfile(catch_path):
                try:
                    #eeeee
                    data_file = load(catch_path,allow_pickle=True,encoding='latin1')
                    tvec = data_file['tvec']
                    data = data_file['data']
                    data_err = data_file['data_error']
                    if 'wrong_det' in data_file:
                        self.wrong_dets_damaged = data_file['wrong_det']
                        if len(self.wrong_dets_damaged) > 0 and not isinstance(self.wrong_dets_damaged[0],str):
                            #convert to utf8
                            self.wrong_dets_damaged = [w.decode('utf-8') for w in self.wrong_dets_damaged]

                
                    self.cache  = tvec,data,data_err
                except Exception as e:
                    print(e)
        
        if hasattr(self,'cache'):
            tvec,data,data_err  = self.cache
            imin,imax = tvec.searchsorted([tmin,tmax])
            t_ind = slice(imin,imax+1)
            return tvec[t_ind],data[t_ind],data_err[t_ind]
        
      

        #number of parallel downloads
        numTasks = 8
        
        #only these signals are stored by the slow diagnostics
        # old shots have all slow channels digitized
        if self.toroidal or self.shot < 160000 or self.fast_data_downsapled: 
            index = list(range(self.nl))

        else: # new shots have only ch 1 - 28 digitized slow in each array
            index = r_[0:28, 32:60]
            self.wrong_dets = r_[28:32, 60:64]
            
        t = time.time()
        #print 'loading data'

        tvec = self.tvec
        nt = len(tvec)
        TDIcalls = []

        if self.fast_data_downsapled:
            TDIcall = '\\SX%sF:SX%sF%.2d'
            for cam, channels in self.detectors_dict.items():
                if not self.toroidal: cam = cam[:-1]
                TDIcalls += [TDIcall%(cam,cam,int(ch.split('_')[-1])) for ch in channels]
            if fetch_errors:
                TDIcalls += ['error_of(%s)'%tdi for tdi in TDIcalls]
            TDIcalls = array(TDIcalls) 
            #print(TDIcalls)
            
        else:
            TDIcall = 'PTDATA2("SX%sS%s",%d)'
            for cam, channels in self.detectors_dict.items():
                if not self.toroidal: cam = cam[:-1]
                TDIcalls += [TDIcall%(cam,ch.split('_')[-1], self.shot) for ch in channels]
            TDIcalls = array(TDIcalls)[index]
            
        
        MDSserver = self.MDSconn.hostspec
        MDStree = 'spectroscopy'

        if self.fast_data_downsapled:
            out = mds_par_load(MDSserver, TDIcalls,  numTasks, MDStree, self.shot)
        else:
            out = mds_par_load(MDSserver, TDIcalls,  numTasks)
            
        print('data loaded in %.2f'%( time.time()-t))
        
        from IPython import embed
 
        data = zeros((nt, self.nl), dtype='single')
        data_error = zeros_like(data)+1e-6
        
        #embed()

        for i,ch in enumerate(index):
            data[:size(out[i]),ch] = out[i]
            data[size(out[i]):,ch] = out[i][-1]
            if self.fast_data_downsapled and not isinstance(out[self.nl+i],str) and fetch_errors:
                data_error[:size(out[i]),ch] = out[self.nl+i]
            data_error[size(out[i]):,ch] = infty  #data from TA system are 2s longer than from PA system
        
        data_error[:] = abs(data)*0.02+data.max(1)[:,None]*0.01
        
        #embed()
        #try:
            #print('Saving raw data RM1 plot', )
            #plot(tvec, data[:,32:])
            #axhline(4.8/10,ls='--')
            #savefig(r'C:\Users\odstrcil\tomography\%d.png'%self.shot)
            #clf()
        #except:
            #pass
                
        from scipy.stats.mstats import mquantiles
        #embed()


        #print(', '.join(['%.2f'%mquantiles(abs(d),.999) for d in data.T]))
        if not self.toroidal:
            overburned = abs(data) > 4.6  #only a rough value, different for each LOS
            ' SLOW CH is saturated pos.( >4.95V)at this time'
            
        else:
            
            overburned = abs(data) > 4.6  #only a rough value, different for each LOS
            overburned[:,-12:] = data[:,-12:]>array((5.5,)*4+(4,)*8)[None,:]

        
        #data[overburned] = nan
        

        if self.fast_data_downsapled:
            offset = data[tvec.searchsorted(0.2):tvec.searchsorted(0.6)].mean(0)
            offset = data[tvec.searchsorted(0.):tvec.searchsorted(0.2)].mean(0)
            
        else:
            offset = data[tvec.searchsorted(tvec[-1]-.1):].mean(0)#BUG can failure if shot is too long

        data-= offset

        self.hardcoded_corrections(tvec, data,data_error)

        
        use_dets = ~in1d(self.dets, self.wrong_dets)
        offset = tvec < .1


        if self.filter_data and not self.fast_data_downsapled:
            #filtering of the data from slow breamch of SXR DAS 
            
            # cameras with the same noise 
        
            if self.toroidal:
                noise_cam_ind = {'TOR': hstack((self.cam_ind['45R1'], self.cam_ind['195R1']))}
            else:
                noise_cam_ind = {'90RP1': hstack((self.cam_ind['90RP1a'], self.cam_ind['90RP1b'])),
                                '90RM1': hstack((self.cam_ind['90RM1a'], self.cam_ind['90RM1b']))}
                noise_cam_ind = self.cam_ind
                noise_cam_ind = {'90RP1a': hstack((self.cam_ind['90RP1a'], self.cam_ind['90RP1b'])),
                                '90RM1': hstack((self.cam_ind['90RM1a'], self.cam_ind['90RM1b']))}



                
            noise_cam_ind = {c:slice(ii[0],ii[-1]+1) for c,ii in noise_cam_ind.items()} 

            data_ = (data[offset]-data[offset].mean(0))
            filtered_data = copy(data)
 
            from scipy import signal
            fnq = (nt-1)/(tvec[-1]-tvec[0])/2

            #BUG correlated noise between 0 and 500Hz survives!!!
            fmax = 50
            b, a = signal.butter(4, fmax/fnq, 'low')
            
            
            from scipy.fftpack import next_fast_len 
            fmax = 500
            G = None #FFT of the deconvolution kernel
      
      
            if fnq > fmax and not  self.fast_data_downsapled:
                nt_ = next_fast_len(nt)
                #filter tunned to match the real low pass filter of the slow digitizer
                bd, ad = signal.ellip(8, 0.01, 100, fmax/fnq*0.97, 'low')
                h = zeros(nt_)
                h[0] = 1
                H = np.fft.rfft(signal.lfilter(bd,ad,h))
                SNR = 100
                #make deconvolution and filtering kernel
                G = conj(H)/(H*conj(H)+1./SNR)
            
            for cam, ind in noise_cam_ind.items():
                if not self.toroidal and False:
                    ind_finit = any(isfinite(data_error),1)
                    err = data_error[ind_finit,ind].mean(0)*(data_[:,ind].std(0)+1e-6)
                    u,s,v = linalg.svd(data_[:,ind].T/err[:,None], full_matrices=False)#SVD veighted be prescribed error
                    u[isfinite(err),0] *= err[isfinite(err)]
                    
                    u2,s2,v2 = linalg.svd(data_[:,ind].T, full_matrices=False)#SVD veighted be prescribed error
      
                    noise =  (data[:,ind][:,-1]/  u2[-1,0])/s2[0]
                    filtered_data[:,ind] -= outer(u2[:,0]*s2[0], noise ).T #outer(noise, u[:,0])

                
                if G is not None:
                    #signal deconvolution and low pass filtering
                    fsig_filt = np.fft.rfft(filtered_data[:,ind],nt_, axis=0)
                    filtered_data[:,ind] = np.fft.irfft(G[:,None]*fsig_filt,axis=0)[:nt]
              

            data = filtered_data 
            
            if not self.fast_data_downsapled:
            #downsample to 1kHz
            
                downsample = max(1,int(round(fnq/500)))
                nt = nt//downsample*downsample
                if downsample >1:
                    tvec = tvec[:nt].reshape(-1,downsample).mean(1)
                    data_error = data_error[:nt].reshape(-1,downsample, data.shape[1]).mean(1)+data[offset].std(0)/2
                    data_error+= (diff(data[:nt].reshape(-1,downsample, data.shape[1]),2,axis=1)).std(1)*3
                    data = data[:nt].reshape(-1,downsample, data.shape[1]).mean(1)
                    overburned = overburned[:nt].reshape(-1,downsample, data.shape[1]).mean(1)>0.1
            
        if not self.fast_data_downsapled:
            #remove slope 
            i1 = tvec < 0
            i2 = tvec > tvec[-1]-.5
            b1,b2 = data[i1].mean(0), data[i2].mean(0)
            a1,a2 = tvec[i1].mean(), tvec[i2].mean()
            

            data -= ((b2-b1)/(a2-a1)*(tvec[:,None]-a1)+b1)


        offset_err = sqrt(mean(data[tvec<.3]**2,0))/sqrt(sqrt(sum(tvec<.3)))
        data_error +=  offset_err[None,:]  #factor 2 is just manul tunning

        cam_ind = {c:slice(ii[0],ii[-1]+1) for c,ii in self.cam_ind.items()} 

        for cam, ind in cam_ind.items():
            
            if cam in self.calib:
                data[:,ind]  *= self.calib[cam][:ind.stop-ind.start]
                data_error[:,ind] *= self.calib[cam][:ind.stop-ind.start]
            else:
                print('Warning: calibration for camera %s is not availible'%cam)
  
        #data_error[overburned] = infty
        
        
        #trick to replace missign fast channel by a slow data 
        if all(std(data[:,[15,20]],0)<100):
            try:
                slow_data = load(catch_path+'_')
                data[:,15] = np.interp(tvec,slow_data['tvec'],slow_data['data'][:,15])
                data[:,20] = np.interp(tvec,slow_data['tvec'],slow_data['data'][:,20])
                data_error[:,15] = np.interp(tvec,slow_data['tvec'],slow_data['data_error'][:,15])
                data_error[:,20] = np.interp(tvec,slow_data['tvec'],slow_data['data_error'][:,20])
            except:
                pass
        
        # data for to accelerate further loading
        self.cache = tvec, data,data_error
        
        ##store data for faster loading 
        if config.useCache:
            savez_compressed(catch_path,tvec=tvec, data=data, data_error=data_error,
                                wrong_det = self.wrong_dets_damaged)

        imin,imax = tvec.searchsorted([tmin,tmax])
        t_ind = slice(imin,imax+1)
        self.tvec = tvec
        
 

        return tvec[t_ind], data[t_ind],data_error[t_ind] 
        

 

        

    def hardcoded_corrections(self,tvec, data,data_err):
        #print 'hardcoded_corrections'
        if not self.toroidal:
            
            
            #data_err[:,self.cam_ind['90RM1a'][2]] = infty
            if self.shot < 166887:  #od 167437 je to zase?
                data[:,self.cam_ind['90RP1a'][15]] *= 1.1

            
            ch1 = self.cam_ind['90RM1a'][9]
            ch2 = self.cam_ind['90RM1a'][10]

            if not (self.fast_data or self.fast_data_downsapled) :
                data[:,[ch1,ch2]]     =  data[:,[ch2,ch1]]
                data_err[:,[ch1,ch2]] =  data_err[:,[ch2,ch1]]
                           
       
            
            
            #signal is inversed? 
            ch = self.cam_ind['90RP1b'][5]
            data[:,ch]*= sign(mean(data[:,ch]))
            data_err[:,ch]*=100
            
       
            
            #NOTE from sxrInvrt2.m from Eric Hollmann
            # channels which need to be inverted for 2007 data
            # dataArray(:,32+9) = -1*dataArray(:,32+9);
            # dataArray(:,32+10) = -1*dataArray(:,32+10);
            # dataArray(:,32+12) = -1*dataArray(:,32+12);
            # dataArray(:,32+16) = -1*dataArray(:,32+16);
            # dataArray(:,32+23) = -1*dataArray(:,32+23);

            # RMP day where RM1 array appears about 2x weaker than RP1 array for some reason
            if (self.shot>154012) and (self.shot<154032):
                data[:,self.cam_ind['90RP1a']]/= 2
                data[:,self.cam_ind['90RP1b']]/= 2

            
            # VDE day where RM1 array still appears about 2x weaker than RP1 array
            if (self.shot>154134) and (self.shot<154161):
                data[:,self.cam_ind['90RP1a']]/= 2
                data[:,self.cam_ind['90RP1b']]/= 2
            
            # RE plateau day where RM1 array also appears about 2x weaker than RP1 array
            if (self.shot>158853) and (self.shot<158877):
                data[:,self.cam_ind['90RP1a']]/= 2
                data[:,self.cam_ind['90RP1b']]/= 2
          
            # all data channels inverted for 2008 data
            if (self.shot<136160) and (self.shot>131000):
                dataArray *=-1;
                
            # upper array with /4 gain at this point
            if (self.shot > 142034) and (self.shot < 142734) : 
                data[:,self.cam_ind['90RP1a']]*= 4
                data[:,self.cam_ind['90RP1b']]*= 4
            
    
        
        if self.toroidal:
            #BUG od 166566 je to uz vetsi korekce  167451 to zase sedi? 
            data[:,self.cam_ind['45R1'][0]] *= 1.17
            data_err[:,self.cam_ind['45R1'][0]] *= 1.17
            
            
            if self.shot < 168847:
                data[:,self.cam_ind['195R1'][0]] *= 2*0.96  #od 169350 je to horsi, zeby uz od 168847, od 169593zase prestrelene
                data_err[:,self.cam_ind['195R1'][0]] *= 2*0.96
            elif self.shot < 175669:# 169593:
                data[:,self.cam_ind['195R1'][0]] *= 2*0.96/1.53
                data_err[:,self.cam_ind['195R1'][0]] *= 2*0.96/1.53
                
            elif self.shot <  181100:
                data[:,self.cam_ind['195R1'][0]] *= 2*0.96/1.53/1.1
                data_err[:,self.cam_ind['195R1'][0]] *= 2*0.96/1.53/1.1
            else:
                data[:,self.cam_ind['195R1'][0]] *= 2*0.96/1.53/1.1*1.1
                data_err[:,self.cam_ind['195R1'][0]] *= 2*0.96/1.53/1.1*1.1
                data[:,self.cam_ind['195R1'][1]] *=   1.1
                data_err[:,self.cam_ind['195R1'][1]] *=   1.1
                
            
            if self.shot > 170000: #just a guess!!      
                data_err[:,self.cam_ind['45R1'][10]] *= 1000 #corrupte channel
            
            if self.fast_data or self.fast_data_downsapled:
                ch = self.cam_ind['195R1'][-2]
                #print self.toroidal, ch, sign(mean(data[:,ch]))
                data[:,ch] *= -1 #BUG? in shot 172221

            #etendue_correction = array([1.251,1.140,1.051,0.985,0.943,
                        #0.928,0.940,0.982,1.054,1.158,1.295,1.467])
            #data[:,self.cam_ind['45R1']]*= etendue_correction
            #data[:,self.cam_ind['195R1']]*= etendue_correction

            if self.shot > 163300 and self.shot < 168847: #just guess
                self.wrong_dets_damaged = ('45R1_11',)
     
            if self.shot > 168847 and not (183200 < self.shot < 183214):      
                self.wrong_dets_damaged = self.detectors_dict['45R1']
                
                     #catch_path = self.geometry_path+'/data/%d_%s.npz'%(self.shot, suffix)
            #new preamplifiers in the 2024 campaign
            if self.shot > 198650: 
	        data[:,self.cam_ind['195R1']] *= -1
		data[:,self.cam_ind['45R1']] *= -1

    
 
def min_fine(x,y): #find extrem with subpixel precision
    i = argmin(y)
    if i==0 or i == len(y)-1: return y[i],x[i]    #not working at the edge! 
    A = x[i-1:i+2]**2, x[i-1:i+2], ones(3)
    a,b,c = linalg.solve(array(A).T,y[i-1:i+2])
    
    ymin =  c-b**2/(4*a)
    xmin = -b/(2*a) 
    return  ymin,xmin


    
     
def main():
    
    
    
    #r_p,z_p,r2_p,z2_p,xangle_p = xraypaths(1e6, toroidal=False)
        
    
    mds_server = "atlas.gat.com"
    #mds_server = "localhost"

    import MDSplus as mds
    MDSconn = mds.Connection(mds_server )
    
    shots = [147008,147010,147013,147014,147015,147016,147017,147018,147019,147022,147023,147024,147025,147026,147027,147028,147029,147030,147032,147033,147034,147035,147036,147043,147044,147045,147046,147047,147048,147051,147055,147059,147062,147063,147066,147067,147068,147069,147070,147071,147072,147078,147079,147080,147082,147084,147092,147094,147095,148775,148776,148777,148778,148780,148781,148783,148784,148785,148786,148787,148788,148789,148790,148791,148792,148793,148794,148796,148797,148798,148799,148800,148801,148802,148803,148804,149392,149393,149394,149395,149396,149397,149398,149399,149400,149401,149403,149404,149405,149406,149407,149408,149409,149410,149411,149412,149413,149414,149415,149416,149682,149683,149685,149686,149687,149688,149689,149691,149692,149693,149694,149695,149696,149697,149698,153132,153133,153134,153135,153136,153137,153139,153140,153141,153142,153143,153144,153145,154971,154976,154979,154980,154981,154982,154983,154984,154985,154986,154987,154988,154989,154990,154991,154992,154993,154994,154995,154996,154997,157444,157445,157446,157447,157448,157449,157450,157451,157452,157453,157454,157455,157456,157457,157458,157459,157460,157461,157462,157531,157813,157814,157815,157817,157818,157819,157820,157836,157838,157841,157858,157859,157860,157862,157865,157866,157867,157870,157871,160093,160094,160095,160096,160099,160101,160102,160103,160104,160123,160147,160149,168969,168970,168971,168972,168973,168974,168976,168977,168980,168981,168982,168983,168984,168985,168986,168987,168988,168989,168990,168991,170478,170479,170480,170481,170483,170484,171222,171223,171224,171226,173608,173609,173610,173611,173612,173613,173614,173615,173618,173621,173622,173623,173624,173625,173627,173628,173630,173632,173633,173634,173635,173636,173637,173642,173643,173646,173647,173648,173649,174480,174481,174482,174483,174486,174487,174488,174490,174491,174492,174493,175302,175305,175309,175313,180447,180448,180455,180458,180461,185211,185212,185213,185306,185312,185314,185332,185333,185334,185344,185352,185353,185356,185399,185403,185404,185408,185412,186600,186601,186602,186603,186604,186605,186606,186607,186608,186609,186610,186611,186612,186613,186614,186939,186953,187082]
    for shot in arange(173600, 193000):
     
        t = MDSconn.get('PTDATA2("SX90PF20_1",%d)'%shot)

        if np.size(t.data()) < 2:
            continue
        MDSconn.openTree('D3D', shot)
        print(shot, MDSconn.get(r'\D3D::TOP.COMMENTS:BRIEF').data())
        
        
        
        
    
    config.useCache=False
    shot = 181536
    x = MDSconn.get('PTDATA2("X1VOLTS",%d)'%shot).data()
    t = MDSconn.get('dim_of(PTDATA2("X1VOLTS",%d))'%shot).data()

    
    PTNAME = '_x=PTDATA2("LYA1%s%.2dRAW",%d,1)'
 
    x2 = MDSconn.get(PTNAME%('H',1,shot)).data()
    t2 = MDSconn.get('dim_of(_x)').data()
  
    x_ = (x-mean(x))/std(x)
    x2_ = (x2-mean(x2))/std(x2)
        
    import IPython
    IPython.embed()
    
    
    #from scipy.signal import lfilter
    
    #a = .05
    #x3 = lfilter([a],[1,-1+a],x_)
    ##plot(t,x3/std(x3))
    ##plot(t, x_)
    
    #plot(t2,x2_ )
    #plot(t2, x2_-interp(t2-0.02,t,x3/std(x3) ))
    ##plot(t, x)
    #xlim(5000,5001)
    #show()
    
    #fx=np.fft.rfft(x_)
    #fx2=np.fft.rfft(x2_)
    
    ind =(t > 0)&(t < 7000)
    x_ = x_[ind]
    
    ind = (t2 >= 0)&(t2 < 7000)
    x2_ = x2_[ind]
    
    
    fx=np.fft.rfft(x_)
    fx2=np.fft.rfft(x2_)
    fx_ = zeros_like(fx2)
    fx_[:len(fx)] = fx
    
    plot(abs(fx)/len(x_)) 
    plot(abs(fx2)/len(x2_)) 
    show()
    
    
    x2r = reshape(x2, (100,-1))
    
    x2r = reshape(x2, (-1,500))

    
    

    #sxr_fast = loader_SXR(shot, '/home/tomas/tomography/geometry/DIIID/SXR/',MDSconn, True, False)

    #tvec_fast, data_fast, _ = sxr_fast.get_data(0,6)
    
    sxr_fast = loader_SXR(shot, '/home/tomas/tomography/geometry/DIIID/SXR/',MDSconn, False, False)

    tvec_slow, data_slow, _ = sxr_fast.get_data(0,6)
    
    import IPython
    IPython.embed()
    
    
    print( MDSconn.get('PTDATA2("SX90PF20_0",%d)'%174720 ))

    #for shot in range(175000, 173902,-1):
        #t = MDSconn.get('PTDfATA2("SX45F01_3",%d)'%shot)
        #if len(t) > 2: print(shot)
        #else: 
            #print 'no ', shot 
    #shot = 175860
    #exit()
    
    shots = loadtxt('/home/tomas/Desktop/DIII-D/carbon_database/database1.txt').T[0]
    #154858
    #sxr_ta = loader_SXR(shot, '/home/tomas/tomography/geometry/DIIID/SXR/',MDSconn, False, True)
    #tvec_ta, sxr_ta, _ = sxr_ta.get_data(0,6)
    try:
        os.mkdir('sxr_database')
    except:
        pass
    for shot in arange(176979, 183000):
         #import IPython
        #IPython.embed()
        t = MDSconn.get('PTDATA2("SX90PF20",%d)'%shot)
        try:
            print(shot, t)
        except:
            continue
        if np.size(t.data()) < 2:
            continue
        
        
        

        sxr_pa = loader_SXR(shot, '/home/tomas/tomography/geometry/DIIID/SXR/',MDSconn, False, False)
        tvec_pa, sxr_pa, _ = sxr_pa.get_data(1,6)
        N = len(tvec_pa)/1000

        imshow(sxr_pa[:len(sxr_pa)/N*N].reshape(-1,N,sxr_pa.shape[1]).mean(1),aspect='auto',interpolation='nearest', vmin=0, origin='lower')
        savefig('sxr_database/%d.png'%shot)
        clf()
        
        
    #sxr_pa = sxr_pa[:,:32]

    exit()
    
    SXR = hstack((sxr_pa[:,:32],sxr_ta[:,:12]))
    
    cov_mat_ = corrcoef(SXR,rowvar=False)
    cov_mat_[32:,:-12]
    
    COV = copy(cov_mat_[:-12,32:])
    COV[10:,2:] = 0
    COV[~isfinite(COV)] = 0
    imshow(COV, vmin=0.9,vmax=1, cmap='seismic', interpolation='nearest', origin='lower',aspect='auto');colorbar();show()
    
    
    x = array([min_fine(arange(32),-c) for c in COV.T])[:,1]
    x_ = polyval(polyfit(arange(12),x,3),arange(12))
    
    plot(SXR.mean(0))
    plot(x_,sxr_ta[:,:12].mean(0) )

    from scipy.interpolate import interp1d
    
    y_ = interp1d(arange(32), SXR[:,:32].mean(0),kind='quadratic')(x_)
    plot(x_,y_)
    
    
    
    eten_correction =  y_/sxr_ta[:,:12].mean(0)
    eten_correction[-1] = 1.5
    eten_correction_ = polyval(polyfit(arange(12),eten_correction,3),arange(12))
    plot(eten_correction)
    plot(eten_correction_)
    show()
    
    

    

    
    #from matplotlib.pylab import *
    ion()
    
    plot(tvec_fast.reshape(-1,25).mean(1),data_fast.reshape(-1,25,data_fast.shape[1]).mean(1)[:,9] )
    
    gca().set_color_cycle(None)

    plot(tvec_slow, data_slow[:,9],'--')
    xlim(5.13, 5.15)
    show()
    
    
    
    
    tvec_ = tvec_fast.reshape(-1,25).mean(1)
    sig_orig = data_fast.reshape(-1,25,data_fast.shape[1]).mean(1)[:,9]
    sig_filt = interp(tvec_, tvec_slow, data_slow[:,9])
    imin,imax = tvec_.searchsorted((4.95,4.98))
    
    from scipy import signal

    fmax = 500
    fnq = (len(tvec_)-1)/(tvec_[-1]-tvec_[0])/2  

    b, a = signal.cheby1(7,0.15, fmax/fnq, 'low')
    #b, a = signal.cheby1(7,0.15, fmax/fnq, 'low')
    b, a = signal.butter(9, fmax/fnq, 'low')
    
    
    
    b, a = signal.ellip(8, 0.01, 100, fmax/fnq*0.97, 'low')


    fsig = signal.lfilter(b,a,sig_orig)
    plot(tvec_[imin:imax], sig_filt[imin:imax])
    plot(tvec_[imin:imax], fsig[imin:imax],'--')
    plot(tvec_[imin:imax], sig_orig[imin:imax])
    xlim(4.95,4.98);show()
    
    

    
    
    
    plot(tvec_,deconv,'--')
    plot(tvec_,sig_filt)
    plot(tvec_,sig_orig)
    ylim(0,25e3)
    xlim(4.95,4.98)
    
    
    
    

    
    #plot(fx)
    plot(np.fft.irfft(np.fft.rfft(sig_filt)/np.fft.rfft(fx)))
        

    X,Y = signal.deconvolve(sig_filt, fx)
            
        
    
    plot(abs(np.fft.rfft(sig_orig[imin:imax])))
    plot(abs(np.fft.rfft(sig_filt[imin:imax])))
    
    
    
    

    ffilt =np.fft.rfft(sig_filt[imin:imax])/np.fft.rfft(sig_orig[imin:imax])
    ffilt[10:] = 0
    plot(np.fft.irfft(ffilt))
    
    #MDSserver = self.MDSconn.hostspec
    


    #if fast_data:
        #print 'fetching time vec'
        #cam = cams[0][:-2] if toroidal else '90'+cams[0][3]
        #ch = 1
        #TDI = ['dim_of(PTDATA2("SX%sF%.2d_%d",%d,1))'%(cam,ch,i,self.shot) for i in range(n_chunks)]
        #self.tvec_fast = mds_par_load(MDSserver, TDI,8)
        #if len(self.tvec_fast[0]) == 1:
            #raise Exception("fast data are not availible!")
        
     
if __name__ == "__main__":
    main()
    
173902
173901
173900
173899
173898
173897
173896
173895
173894
173893
173892
173891
173890
173889
173888
173887
173886
173885
173884
173883
173882
173881
173880
173879
173878
173877
173876
173875
173874
173873
173872
173871
173870
173869
173868
173867
173866
173828
173827
173826
173825
173824
173823
173822
173821
173820
173819
173818
173817
173816
173815
173814
173813
173812
173811
173810
173809
173808
173807
173806
173805
173804
173803
173802
173801
173800
173799
173798
173683
173682
173681
173680
173679
173678
173677
173676
173675
173674
173673
173672
173671
173670
173669
173668
173667
173666
173665
173664
173663
173662
173661
173660
173659
173658
173657
173656
173655
173654
173650
173649
173648
173647
173646
173645
173644
173643
173642
173641
173640
173639
173638
173637
173636
173635
173634
173633
173632
173631
173630
173629
173628
173627
173626
173625
173624
173623
173622
173621
173620
173619
173618
173617
173616
173615
173614
173613
173612
173611
173610
173609
173345
173344
173343
173342
173341
173340
173339
173338
173337
173336
173335
173334
173333
173332
173331
173330
173329
173328
173327
173326
173323
173322
173321
173320
173319
173318
173317
173316
173315
173314
173313
173312
173311
173310
172721
172720
172718
172351
172350
172349
172348
172347
172346
172345
172344
172343
172342
172341
172340
172339
172338
172337
172336
172335
172334
172333
172332
172331
172330
172329
172328
172327
172326
172325
172323
172322
172321
172320
172315
172314
172313
172312
172311
172310
172309
172308
172307
172306
172305
172304
172303
172302
172301
172300
172299
172298
172297
172296
172295
172294
172293
172292
172286
172285
172284
172283
172282
172281
172280
172279
172278
172277
172276
172275
172274
172273
172272
172269
172268
172267
172266
172265
172260
172259
172258
172257
172256
172255
172254
172253
172252
172251
172250
172249
172248
172247
172246
172245
172244
172243
172242
172241
172240
172239
172238
172237
172229
172228
172227
172226
172225
172224
172223
172222
172221
172220
172219
172218
172217
172216
172215
172214
172209
172208
172207
172206
172205
172204
172203
172202
172201
172200
172108
172107
172106
172105
172104
172103
172102
172101
172100
172099
172098
172097
172096
172095
172094
172093
172092
172069
172068
172067
172066
172065
172064
172063
172062
172061
172060
172059
172058
171616
171615
171614
171613
171612
171611
171610
171609
171608
171607
171606
171605
171604
171603
171602
171601
171600
171599
171598
171597
171596
171595
171594
171593
171592
171591
171590
171589
171588
171546
171545
171544
171543
171542
171541
171540
171366
171365
171364
171363
171362
171361
171360
171359
171358
171357
171356
171355
171354
171353
171352
171351
171350
171349
171348
171347
171346
171345
171344
171343
171342
171341
171340
171339
171338
171337
171336
171335
171334
171333
171094
171093
171092
171091
171090
171089
171088
171087
171086
171085
171084
171083
171082
171081
171080
171079
171078
171077
171076
171075
171074
171073
171072
171071
171070
171069
171068
171067
171066
171065
171064
170710
170709
170708
170707
170706
170705
170704
170703
170702
170701
170700
170699
170698
170697
170696
170695
170694
170693
170692
170691
170544
170543
170542
170541
170540
170539
170538
170537
170536
170535
170534
170533
170532
170531
170530
170529
170528
170527
170526
170525
170524
170523
170522
170521
170520
170519
170518
170517
170515
170514
170513
170512
170511
170510
170509
170508
170507
170506
170505
170504
170503
170502
170501
170500
170499
170498
170497
170261
170260
170259
170258
170257
170256
170255
170254
170253
170252
170251
170250
170249
170248
170247
170246
170245
170244
170243
170242
170241
170240
170239
170238
170237
170236
170235
170234
170233
170232
170137
170136
170135
170134
170133
170132
170131
170130
170129
170128
170127
170126
170125
170124
170123
170122
170121
170120
170119
170118
170117
170116
170115
