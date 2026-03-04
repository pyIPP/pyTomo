 #!/usr/bin/env python 
# -*- coding: utf-8 -*-

import sys,os
from time import time
from scipy.interpolate import UnivariateSpline
from matplotlib.pylab import *
from numpy import *

c = 'r','g','b','y','k','m'
s = '*','x','o','+'
files = os.listdir("./")
files.sort()
for i,file in enumerate(files):
    if file.startswith("geom_corrections_"):
        print(file[-9:-4])
        try:shot = int(file[-9:-4])
        except: continue
        #i = random.randint(1009)
        cams,calb = loadtxt(file, dtype={'names': ('cam', 'calib'),
                                'formats': ('S4', 'd')}, unpack=True)
        plot(calb,s[i%len(s)],color=c[i%len(c)],label=shot)
        
legend()
show()

