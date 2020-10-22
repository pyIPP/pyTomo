#!/usr/bin/env python 
# -*- coding: utf-8 -*-

import sys,os
from time import time
from scipy.interpolate import UnivariateSpline

from matplotlib.pylab import *
from numpy import *



#==================  Parameters  =========================================

shot = 34620
tmin = 5.6 
tmax = 6.5


#shot = 34695
#tmin = 1.8 
#tmax = 4
shot = 36476
tmin = 1
tmax = 6




tokamak = 'ASDEX'
diagnostic = 'SXR'
original_corrections = '/afs/ipp/home/t/todstrci/pytomo_orig/geometry/%s/%s/det_pos_corr_2019_guess.txt'%(tokamak,diagnostic)
modified_corrections = '/afs/ipp/home/t/todstrci/pytomo_orig/geometry/%s/%s/det_pos_corr_2019_1.txt'%(tokamak,diagnostic)
dets = ['G','F','H1','H2','H3','I2','J1','J2','J3','K1','K2','L','M']
#dets = ['G','F','H2','I2','J2','K1','K2','L','M']


#==================  tomography setting  =========================================

command = 'pytomo --lambda_solver .7 --tmin %.3f --tmax %.3f -B 4 -x 80 -y 120  -d 100 -u 100  -T 3 --transform_order_a 1 --transform_order_r 50 -S 0 --ratiosolver 3 --solver 3 %d '
#command = 'pytomo --lambda_solver .7 --tmin %.3f --tmax %.3f -B 16  -d 100 -u 100   --ratiosolver 3  %d '


#python ./pytomo.py 34695 --lambda_solver 0.55 --tmin 1.8  --tmax  4  -B 4 -x 80 -y 120  -d 100 -u 100  -T 3 --transform_order_a 0 --transform_order_r 50 -S 0 --ratiosolver 3 --solver 3 

#python ./pytomo.py 34620 --lambda_solver 0.55 --tmin 5.6  --tmax  6.5  -B 4 -x 80 -y 120  -d 100 -u 100  -T 3 --transform_order_a 0 --transform_order_r 50 -S 0 --ratiosolver 3 --solver 3 

pos_corr = loadtxt(original_corrections,
            dtype={'names': ('det', 'alpha'),'formats': ('U2',  'f4')})
pos_corr =  {k:float(item) for k,item in pos_corr}



from os.path import expanduser 

#fig = figure(shot,figsize=(12, 8))
#fig.subplots_adjust(hspace=0.15, wspace = 0.1)
#fig.suptitle('Optimization of the detectors position at shot %d'%shot)
#fig.show()


def find_chi2(calib_dict):
    f = open(modified_corrections, 'w')

    for k,i in list(pos_corr.items()):
        f.write('%s  %.4f\n'%( k,i ))
    f.close()

    ##print(command%(tmin,tmax,shot))
    os.system(command%(tmin,tmax,shot))
    #exit()

    

    chi2 = load(expanduser('~/tomography/tmp/convergence_%d.npz'%shot))['chi2']
    #print('chi2', chi2)
    #chi2 = sort(chi2)[:-len(chi2)//10]  #remove 10% of the most extreme values
    
    #print(exp(mean(log(chi2))))
    #exit()
    return exp(mean(log(chi2)))
    

cams,calib = loadtxt('/afs/ipp/home/t/todstrci/tomography/geometry/%s/%s/calibration/%d.txt'%(tokamak,diagnostic, shot), 
                     dtype={'names': ('cam', 'calib'),'formats': ('U4', 'd')}, unpack=True)





#nrow = max(len(dets)/4,1)
#ncol = int(ceil(len(dets)/float(nrow)))


for jd, j in enumerate(dets):
    print(j, ' %.3f'%pos_corr[j])
    
    
    
def get_grad(pos_corr,da):
    
    chi2_0 = find_chi2(pos_corr)
    Grad = zeros(len(dets))
    

    for jd, j in enumerate(dets):
        pos_corr[j] += da
        Grad[jd] = (find_chi2(pos_corr)-chi2_0)/da
        pos_corr[j] -= da
        print('gradient', jd, j)
    return Grad
    
def line_search(Grad,pos_corr,da):
    #stupid but robust line optimization
    chi2_0 = find_chi2(pos_corr)
    Grad/= linalg.norm(Grad)
    
    while abs(da) > 0.05:
        for jd, j in enumerate(dets):
            pos_corr[j]-= Grad[jd]*da
        chi2 = find_chi2(pos_corr)
        
        print('============================================= da:%.3f   chi2_init: %.3f    chi2min: %.3f    chi2: %.3f  ,169.7 ==================================='%( da, chi2_init,chi2_0,chi2,))


        if chi2 > chi2_0:
            for jd, j in enumerate(dets):
                pos_corr[j]+= Grad[jd]*da
            da/= -2
        else:
            da*= 1.5
            chi2_0 = chi2
        
    return pos_corr
    


chi2_init = find_chi2(pos_corr)
Grad = get_grad(pos_corr,0.2)
pos_corr = line_search(Grad, pos_corr, 0.2)
Grad = get_grad(pos_corr,0.1)
pos_corr = line_search(Grad, pos_corr, 0.2)
Grad = get_grad(pos_corr,0.05)
pos_corr = line_search(Grad, pos_corr, 0.2)

Grad = get_grad(pos_corr,0.04)
pos_corr = line_search(Grad, pos_corr, 0.2)



print('done')
print('Check the results please')
