import matplotlib.pylab as plt
import numpy as np
import MDSplus


class READ_MDS:


    def __init__(self):

        server = 'localhost:8000'
        #slogin -l todstrci  -L 8000:mdsplus:8000 gate2.aug.ipp.mpg.de 
        #TODO run ssh tunneling ? check if is it in the intranet? 

        #server = 'mdsplus.aug.ipp.mpg.de'
        #conn = MDSplus.Connection(server)
        
        conn = MDSplus.Connection( server)

        print('connected to %s' %tok)

    def cShotNr(self, exper, diag,shot):
        #BUG not finished
        return -1
    #def 
            
    def GetParameter(self, nshot, diag, parameter_set,parameter,exper='AUGD',ed=0):
        #https://www.aug.ipp.mpg.de/aug/local/mdsplus/tdi/augparam.html
        
        mds_str = '_s=augparam (%d,"%s","%s","%s","%s",%d)'%(nshot,diag,parset,par,exper, ed)

        out = conn.get(mds_str)
        if isinstance(out,MDSplus.mdsscalar.String):
            return str(out)
        else:
            return np.asarray(out)

        
    def GetVessel(self, nshot = -1):
 
        mds_str = 'augvessel (%d, _R, _z, _index, _oshot)'%nshot
        conn.get(mds_str)
        R = np.asarray(conn.get('_R'))
        z = np.asarray(conn.get('_z'))
        ind = np.asarray(conn.get('_index'))
        n_poly = len(ind)-1
        structs = []
        for i in range(n_poly):
            structs.append((R[ind[i]:ind[i+1]], z[ind[i]:ind[i+1]]))
            
        #for sx,sy in structs: plot(sx,sy)
            
        return structs
        

    def GetLastShot(self):
        
        mds_str = ' _shot = auglastshot()'
        return int(conn.get(mds_str))
    
    def Eqmapping(self):
        #https://www.aug.ipp.mpg.de/aug/local/mdsplus/tdi/augconv.html
        #not implemented yet!
        


        
        
    def GetSignal(self, nshot, diag, sig,exper='AUGD', ed=0,t1=0,t2=10):
        #https://www.aug.ipp.mpg.de/aug/local/mdsplus/tdi/augdiag.html
 

        mds_str = '_s=augdiag(%d,"%s","%s","%s",%d,%f,%f,_oshot,_oedition)' %(nshot, diag, sig, exper, ed, t1, t2)

        try:
            self.data = np.array(conn.get(mds_str)) #.T
            ndims = self.data.ndim
            self.units = str(conn.get('units_of(_s)'))
            self.shot = int(conn.get('_oshot'))
            self.ed = int(conn.get('_oedition'))

            if ndims == 1:
                self.time = np.array(conn.get('dim_of(_s)'))
            elif ndims == 2:
                self.time = np.array(conn.get('dim_of(_s,0)'))
                self.R    = np.array(conn.get('dim_of(_s,1)'))
            elif ndims == 3: #can be wrong?
                self.R    = np.array(conn.get('dim_of(_s,0)'))
                self.z    = np.array(conn.get('dim_of(_s,1)'))
                self.time = np.array(conn.get('dim_of(_s,2)'))
        except:
            self.error = 1
            
        return self.data 



if  __name__ == "__main__":

#    nshot = 66128
#    dda = 'CXFM'
#    datatype = 'TI' # 'ANGF'
#    mds = READ_MDS(nshot, dda, datatype, tok='JET')

#    nshot = 28053
    diag  = 'CEZ'
    sig   = 'Ti_c'
    nshot = 30301
   diag = 'EQI'
    #diag = 'EQH'
    sig = 'PFM'
    mds = READ_MDS(nshot, diag, sig)

    print 'Time array:', mds.time
    plot = False
    if plot:
        jr = 0
        jt = 10
        plt.subplot(2, 1, 1)
        plt.plot(mds.time, mds.data[:, jr])
        plt.title('R = %8.4f' %mds.R[jr])
        plt.subplot(2, 1, 2)
        ind = (mds.data[jt, :] > 0)
        print ind
        plt.plot(mds.R[ind], mds.data[jt, ind])
        plt.title('t = %8.4f' %mds.time[jt])

        plt.show()
