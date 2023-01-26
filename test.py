#!/usr/bin/env python
import os, sys, time
from multiprocessing import Pool
import numpy as np

z = []    # integers are immutable, let's try mutable object

#global x
###x = {'x':np.zeros(1000000)} # replace `23` due to small integers share representation

#x = {'x':np.ones(1000000)} # replace `23` due to small integers share representation
    
cursor = None
def initializer(xval):
    global x
    x = xval
    
def run():

    
    pool = Pool(processes=4, initializer=initializer, initargs= ({'x':np.ones(1000000)},))
    pool.map(printx, (1,2,3,4))
    
    
    
    
    
def printx(y):
    global x
    #print(x)
    if y == 3:
       x['x'] = -x['x']
    z.append(y)
    print( os.getpid(), sum(x['x']), id(x['x']), z, id(z) )
    print( y)
    if len(sys.argv) == 2 and sys.argv[1] == "sleep":
       time.sleep(.1) # should make more apparant the effect




if __name__ == '__main__':
    run()


    
    
    
    
