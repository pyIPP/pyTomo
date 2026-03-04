import sys, os
import MDSplus as mds
#sys.path.append("/fusion/projectEs/codes/pytomo/geometry/DIIID")
from geometry.DIIID.mag_equ import  Equlibrium
from numpy import *
from geometry.DIIID import map_equ

def main():
    
    shot = int(sys.argv[1])
    mag_diag = sys.argv[2]

    MDSconn = mds.Connection('atlas' )

    eqm = map_equ.equ_map(MDSconn)   
    mag_exp = 'DIIID'
    mag_ed = 0
 

    eq_diags = ((mag_diag,mag_exp,mag_ed),('EFIT01','DIIID',0),('EFITRT1','DIIID',0))

    stat = False

    for diag,mag_exp,mag_ed  in eq_diags:
        stat = eqm.Open(shot, diag=diag, exp=mag_exp, ed=mag_ed)
        if stat and size(eqm.t_eq) >2 : break
        print('Warning: equlibrium for shot:%d diag:%s  exp:%s  ed:%d  was not found!! other will be used'%(shot,diag,mag_exp,mag_ed))

    if not stat:
        raise Exception('equlibrium for shot:%d diag:%s  exp:%s  ed:%d  was not found!! use another one')


    mag_diag = diag

    EQU = Equlibrium(MDSconn,eqm, shot, mag_diag,mag_exp,mag_ed)

    if 'TRA' in mag_diag:
        output = EQU.getTranspEquilibrium()
    else:
        output = EQU.getStandartEquilibrium()
    output[mag_diag] = mag_diag


    tsurf = output['tsurf']
    surf_coeff = output['surf_coeff']
    mag_axis = output
    mag_axis['tvec'] = output['tvec_fast']

    try:
            
        pfm = copy(eqm.pfm)
        pfm-= eqm.psi0
        pfm/= (eqm.psix-eqm.psi0)
        pfm_tvec = eqm.t_eq
        pfm_R = eqm.Rmesh
        pfm_Z = eqm.Zmesh
        PFM = {'pfm':single(pfm),'pfm_tvec':pfm_tvec,'pfm_R':pfm_R,'pfm_Z':pfm_Z}
        output.update(PFM)
    except:
        pass
    
    path = '/local-scratch/'+os.environ.get('USER')+'/'
    try:
        os.mkdir(path)
    except:
        pass
    savez_compressed(path+'/MagField_fast_%d.npz'%shot,diag=mag_diag,exp=mag_exp,ed=mag_ed,**output)


    MDSconn.closeAllTrees()









if __name__ == "__main__":
    main()


