#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is main program file. It calls all other submodules. This file contains default setting
for tomography and command-line interface. All functions should be fully accesible via the commandline.



    :Usage: `start.py [ options ]  [ action ] [ shot ]`
    :Example:
            >>> start.py -r -d 5 -w 1,2,3  --reconstruct 2721
            >>> start.py -v -r -d 5 -w 1,2,3 -S4 -R 2721
            >>> start.py --verbose --rapid --data_smooth=5 --wrong_dets=1,2,3 --regularization=4 -G 2721

            ``start.py  or start.py [ options ] -G [shot]`` => start GUI
    :Options: ``-[adeiGDhnopPRrsvwxytb]``

            ``[ --help  ] [ --regularization=TYPE ] [ --danis=RATIO    ] [ --data_smooth=NBINS]``

            #``[ --debug ] [ --input=INPUT    ] [ --noboundary ] [ --output=OUTPUT    ] [ --sparallel ]``

            ``[ --nx=XPIX ] [ --ny=YPIX ] [ --rapid ] [ --separated-plots  ] [ --errorscale=SCALE ]``

            ``[ --wrong_dets=LIST ] [ --tmin=TMIN   ] [ --tmax=TMAX ] [ --verbose ] [ --enable_output ]``

            ``[ --postprocessing ]  [ --tokamak=NUM ] [ --rapid_blocks= ] [ --plot_svd ] [ --wavelets ]``

    :Actions:  ``[ --preview  | --reconstruct | --gui | --web ]``

    :Smoothing: 1 - Isotropic minimum fisher information regularization

            2 - Anisotropic minimum fisher information regularization

            3 - Isotropic minimum diffusion regularization

            4 - Anisotropic minimum diffusion regularization
            
            4    Energy 
            
            5    Isotropic MNGR
            
            6    Anisotropic MNGR
            
            7    Anisotropic Imgesson




"""








#try:
import matplotlib   #on JET must  be used Qt4Agg backend
matplotlib.rcParams['backend'] = 'Qt4Agg'  #choose one of GTK GTKAgg GTKCairo CocoaAgg FltkAgg MacOSX QtAgg Qt4Agg TkAgg WX WXAgg Agg Cairo GDK PS PDF SVG
##!! with no X use Agg !! 
matplotlib.rcParams['backend'] = 'Agg'  #choose one of GTK GTKAgg GTKCairo CocoaAgg FltkAgg MacOSX QtAgg Qt4Agg TkAgg WX WXAgg Agg Cairo GDK PS PDF SVG

if matplotlib.compare_versions(matplotlib.__version__, '1.9.9'):
# http://matplotlib.org/users/dflt_style_changes.html
    params = { 
            'axes.labelsize': 'medium',
            'axes.titlesize': 'medium',
            'xtick.labelsize' :'medium',
            'ytick.labelsize': 'medium',
            'font.size':12,
            'mathtext.fontset': 'cm',
            'mathtext.rm': 'serif',
            'grid.color': 'k',
            'grid.linestyle': ':',
            'grid.linewidth': 0.5,
            'lines.linewidth'   : 1.0,
            'lines.dashed_pattern' : (6, 6),
            'lines.dashdot_pattern' : (3, 5, 1, 5),
            'lines.dotted_pattern' : (1, 3),
            'lines.scale_dashes': False,
            'errorbar.capsize':3,
            'mathtext.fontset': 'cm',
            'mathtext.rm' : 'serif',
            'legend.loc':'upper right',
            'legend.fontsize':'large',
            'legend.framealpha':None,
            'legend.scatterpoints':3,
            'legend.edgecolor':'inherit'}
                

    
    
    matplotlib.rcParams.update(params)


import sys,os
from numpy import *
from scipy import sparse
from numpy import linalg
import time
from scipy.interpolate import RectBivariateSpline
from scipy.io import loadmat
import pickle
import os
import os.path
import socket
from main import *
from shared_modules import read_config
import config
from shutil import copyfile

set_printoptions(precision=4,linewidth=400)

global cfg_file
cfg_file = "tomography"


def loadSetting( Reload = False):
    """ Subroutine to load saved settings from previous run or default settings if no saved settings exist.

    :param str cfg_file+".npz":  Path to saved setting, generated automatically
    :param bool Reload:         If ``True`` reload all values, else try to load if possible

    :var int shot:      Number of recontructed shot / pulse


    :var bool solid_parallel: Use faster parallel solve, but unbreakable and without progressbar in GUI
    :var int nx:        Number of pixels in X axis
    :var int ny:        Number of pixels in Y axis
    :var int tok_index: Index of choosen tokamak, possible indexes are following:

        
                        ``0-1    Golem`` - bolometry, camera 

                        ``2-4    COMPASS`` - BOLO, SXR, camera

                        ``5-8    JET`` - SXR-slow, SXR-fast, bolo, neutrons

                        ``9-12    ASDEX`` - SXR,SXRfast, bolo,axuv

                        ``13-14    ToreSupra`` - SXR, Camera
                        
                        ``15-20    TCV`` - XTOMO, AXUV, BOLO, DMPX , FIR,XTEPRO
                        ``21-25    DIII-D`` - SXR,  SXR fast,BOLO ,DISRAD

                        

    :var int regularization: Apriory regularization used for minimalization, indexes are following:

                        ``0    MFI``
s
                        ``1    Anisotropic MFI``

                        ``2    MDIFF``

                        ``3    Anis. MDIFF``

                        ``4    Energy ``
                        
                        ``5    Isotropic MNGR``
                        
                        ``6    Anisotropic MNGR``
                        
                        ``7    Anisotropic Imgesson``


    :var bool boundary:    Allow regularization with boundary => The reconstruction will be zero on the boundary
    :var bool plot_all:  Plot each timeslice separately. It can be slow
    :var str local_path: Basic path of the program, default setting is path of start.py
    :var str input_path: Path were are searched input data/geometry for tokamaks
    :var list wrong_dets: List of wrong detectors that are ignored, starts from 0 !
    :var double tmin:     Beginning of reconstruction
    :var double tmax:     End of reconstruction
    :var int data_smooth:  Use moving average for input data over `n` bins
    :var int data_undersampling:    Use every `n`-th snapshot
    :var double danis:    Ratio of anisotropic matrices, B = sigmoid(-n)*Bperp+sigmoid(n)*Bpar  => The bigger number the more will reconstruction follow the field For MFI is recommended ~ 4
    :var int virt_chord: Number of virtual chords used in geometry matrix
    :var double boundary_prec:    Set the pressure forcing reconstruction to be zero on boundary. Precision of the virtual senzor (boundary) is min(sigma)/X. Recommended value is 50 but if it is too unsmooth use lower
    :var double error_scale: Scale the expected errors, The higher number the smoother reconstruction. Original value is 1
    :var bool enable_output: Create plots in printable quality (pdf)
    :var str output_path: Path were are printable plots saved
    :var str hostName:  Name of computer, default is `!uname -n``, is used to remove default setting if is changed
    :var bool gnuplot:  Gnuplot is faster but the output graphs quality are worse
    :var bool SVD:  Substract first mode of SVD from plotted results
    :var bool Chords:   Plot chords to the final recontruction graphs
    :var int ifishmax: Number of Minimum Fisher loops (outer loops), recommended is 3
    :var array magfield_shift:   Artificial shift of magnetic field
    :var bool postprocessing:   Perform postprocessing -- total emissivity, position of center of emissivity, entropy, smoothness
    :var bool gui:       Start graphical interface
    :var bool reconstruct:  Perform tomographic reconstruction
    :var int solver: Used main solver for tomography, 0 - MFI, 1 - Rapid MFI, 2 - Fast SVD,3 -Faster SVD, 4 - QR, 5 - GEV, 6 - GSVD, 7 - None
    :var int presolver: Used main presolver for separated MFI solver, 0 - None, 1 - Rapid MFI, 2 - Fast SVD, 3 - QR, 4 - GEV, 5 - small resolution
    :var int ratiosolver: Used abel like solver for determination of ratio of dets, 0 - None, 1 - Rapid MFI, 2 - Fast SVD, 3 - Faster SVD, 4 - Faster SVD, 5 - QR, 6 - GEV,  7 - small resolution
    :var bool rem_back: Remove "background" -- the first timeslices of the reconstruction
    :var bool positive_constrain:  force solution to be positive (not compatible with solver 0 and 1)
    """


    #os.path.abspath(
    program_path = os.path.abspath(__file__)
    program_path = program_path[:program_path.rfind('/')]+'/'
        
    local_path = os.path.expanduser('~/tomography/')
    if not os.path.exists(local_path):
        try:
            os.mkdir(local_path)
        except:
            print('error: os.mkdir(local_path)')
            raise
    if not os.path.exists(local_path+'/geometry/'):
        try:
            os.mkdir(local_path+'/geometry')
        except:
            print('error: os.mkdir(geometry)')
            raise
        
    if not os.path.isfile(local_path+cfg_file+'.cfg'):
        print('copy original setting')
        copyfile(program_path+cfg_file+".cfg",local_path+cfg_file+'.cfg')
        
    
    inputs = read_config(local_path+cfg_file+".cfg")
 
    inputs['program_path']= program_path
    inputs['local_path']  = local_path
    inputs['output_path'] = os.path.expanduser(os.path.expandvars(inputs['output_path']))+'/'
    inputs['tmp_folder']  = os.path.expanduser(os.path.expandvars(inputs['tmp_folder']))+'/'

        
    config.wrong_dets_pref = inputs['wrong_dets']
    if 'usecache' in inputs:
        config.useCache = inputs['usecache']

    return inputs



help = """
Usage:
        pytomo [ options ]  [ action ] [ shot ]
        Example:
        pytomo  -d 5 -w 1,2,3  --reconstruct=1 30579
        pytomo -v  -d 5 -w 1:5  -S 4 -R 1 -G  30579
        pytomo --verbose  --data_smooth=5 --wrong_dets=1,2,3,4 --regularization=0 -G 30579
        pytomo  or pytomo [ options ] -G [shot] => start GUI
Options: -[adeiGDhnopPRrsvwxytb]
          [ --help  ] [ --regularization=TYPE ] [ --danis=RATIO    ] [ --data_smooth=NBINS]
          [ --debug ] [ --input=INPUT    ] [ --noboundary ] [ --output=OUTPUT    ] [ --sparallel ]
          [ --nx=XPIX ] [ --ny=YPIX ] [ --rapid ] [ --separated-plots  ] [ --errorscale=SCALE ]
          [ --wrong_dets=LIST ] [ --tmin=TMIN   ] [ --tmax=TMAX ] [ --verbose ] [ --enable_output ]
          [ --postprocessing ]  [ --tokamak=NUM ] [ --rapid_blocks= ] [ --SVD ] [--plot_poloidal_spectrum]
Actions:  [ --preview  | --reconstruct | --gui | --web ]

Regularzation:         0 - Isotropic minimum fisher information regularization
                1 - Anisotropic minimum fisher information regularization
                2 - Isotropic minimum diffusion regularization
                3 - Anisotropic minimum diffusion regularization

If an option is not defined, program will use inputs from previous run or default options.
For more information read the source code or contact me at tomas.odstrcil@ipp.mpg.de.
"""

import getopt, sys, os
def usage():

    if sys.platform != 'linux2':
        text.replace('start.py', 'start.bat')
    print(help)
    sys.exit(1)



def main():
    """
    Main function, process options from command-line and calls other submodules.

    :param:  None
    :raises: None

    """
    import argparse
    import matplotlib   #on JET must  be used Qt4Agg backend
    matplotlib.rcParams['backend'] = 'Agg'  #choose one of GTK GTKAgg GTKCairo CocoaAgg FltkAgg MacOSX QtAgg Qt4Agg TkAgg WX WXAgg Agg Cairo GDK PS PDF SVG
    # !! with no X use Agg !! 

 
    inputs = loadSetting( )

    inputs['tmp_folder']  = os.path.expanduser(os.path.expandvars(inputs['tmp_folder' ]))
    inputs['output_path'] = os.path.expanduser(os.path.expandvars(inputs['output_path']))
   
    for f in [inputs['tmp_folder'],inputs['output_path']]:
        if not os.path.isdir( f):
            os.mkdir(f)
            
    if len(sys.argv) == 1:
        inputs['solid_parallel'] = False
        inputs['prune_dets'] = 0
        tok = loaddata(inputs)
        startGUI(inputs, tok)
        sys.exit(0)
        
    elif  len(sys.argv) == 2 and  sys.argv[1] in ('-v', "--verbose", "--debug"):
        config.DEBUG = True

        inputs['solid_parallel'] = False
        inputs['prune_dets'] = 0
        tok = loaddata(inputs)
        startGUI(inputs, tok)
        sys.exit(0)

    

 
    
  
    # default options !!! 
    deriv = False
    plot_projection = False
    
    if 'dpi' not in inputs: inputs['dpi'] = 85
    if 'img_size' not in inputs: inputs['img_size'] = 6


    parser = argparse.ArgumentParser(description='Perform tomographic recontruction', 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,  epilog='''
        Example:
        pytomo  -d 5 -w 1,2,3  --reconstruct=1 30579
        pytomo -v  -d 5 -w 1:5  -S 4 -R 1 -G  30579
        pytomo --verbose  --data_smooth=5 --wrong_dets=1,2,3,4 --regularization=0 -G 30579
        pytomo  or pytomo [ options ] -G [shot] => start GUI

    If an option is not defined, program will use inputs from previous run or default options.
    For more information read the source code or contact me at tomas.odstrcil@ipp.mpg.de''')


    parser.add_argument('shot', metavar='PULSE', help='Number of pulse', type=int, default=inputs['shot'])
    parser.add_argument("-v", "--verbose", "--debug", help=" Set up verbose output on command-line",
                        action='store_true', dest='DEBUG')
    parser.add_argument("-o", "--output", dest='output_path', default=inputs['output_path'],
                        type=str,  help="Path were are printable plots saved")
    parser.add_argument('--enable_output',action='store_true', dest='enable_output',
                        default=inputs['enable_output'], help="Create plots in printable quality (pdf)")
    parser.add_argument("-n", "--noboundary", action='store_const',const=-1,   default=inputs['boundary'],
                        dest='boundary' , help="Do not zero emissivity out of boundary" )
    parser.add_argument("-b", "--boundary",  type=float, default=inputs['boundary'], dest='boundary',
                        help="Boundary regularization (0 = off)" )
    parser.add_argument("-s", "--separated-plots", action='store_true', default=inputs['plot_all'],
                        dest='plot_all', help="Plot each timeslice separately. It can be slow"  )
    parser.add_argument("--solver", dest='solver', default=inputs['solver'] , type=int,
                        help='Used main solver for tomography, \n0 - MFI\n 1 - Rapid MFI\n 2 - SVD,'+\
                        ' 3 - Fast SVD\n 4 - QR\n 5 - GEV\n 6 - GSVD\n 7 - None')
    parser.add_argument("--ratiosolver", dest='ratiosolver', default=inputs['ratiosolver'] ,
                        help='Used abel like solver for determination of ratio of '+\
                            'dets, 0 - None, 1 - Rapid MFI, 2 - SVD, 3 - Fast SVD,'+\
                                '  4 - QR, 5 - GEV,\n 6 - GSVD\n  7 - small resolution', type=int)
    parser.add_argument("--presolver", dest='presolver', default=inputs['presolver'] ,  type=int,
                        help='Used main presolver for separated MFI solver, 0 - None, '+\
                            '1 - Rapid MFI, 2 - SVD, 3 - Fast SVD, 4 - QR, 5 - GEV,'+\
                                '\n 6 - GSVD\n 7 - small resolution')
    parser.add_argument("--lambda_solver", dest='lambda_solver', default=inputs['lambda_solver'] ,
                        type=str, help='Method to find optimal regularization parameter,'+\
                            ' gcv, press, chi2, or any number between 0 and 1 for manual regularization')
    parser.add_argument("--phantom", dest='use_phantom_data', default=inputs['use_phantom_data'] ,
                        type=str, help='Use artificial phantom instead of a real data:'+\
                        ' Asymmetric, Hollow,Gaussian, Hamming,Flat,Peaked ,shepp_logan,Asymmetric_full,snake,mode')
    parser.add_argument("-S", "--regularization", default=inputs['regularization'],  type=int, 
                        help="Apriory regularization used for minimalization, indexes "+\
                        "are following: 0    MFI, 1 Anisotropic MFI, 2    MDIFF, "+\
                            "3    Anis. MDIFF, 4  No derivation, 5 Isotropic MNGR,"+\
                                " 6 Anisotropic MNGR, 7 Anisotropic Imgesson ")
    parser.add_argument("-t", "--tokamak", dest='tok_index',default=inputs['tok_index'],  type=int, 
                help="Index of choosen tokamak, possible indexes are following:0:"+\
                    "Golem - BOLO    1:Golem - Camera    2:COMPASS - BOLO    "+\
                        "3:COMPASS - SXR    4:COMPASS Camera    5:JET - slow SXR "+\
                    "   6:JET - fast SXR  7:JET - neutrons    8:JET - BOLO  "+\
                    "  9:ASDEX SXR    10:ASDEX SXR_fast    11:ASDEX BOLO  "+\
                    "  12:ASDEX AXUV     13:Tore Supra - SXR    14:ToreSupra Camera "\
                    "15 TCV:XTOMO     16 TCV:AXUV     17  TCV:BOLO   18 TCV:DMPX    "\
                    "19 TCV:FIR    20 TCV:XTEPRO    21 DIII-D:SXR   "\
                    "22 DIII-D:SXR fast  23:DIII-D:BOLO  24:DIII-D:DISRAD ")
    parser.add_argument("-B", "--rapid_blocks", type=int, dest='rapid_blocks',
                        default=inputs['rapid_blocks'] , help="number of blocks for rapid methods")
    parser.add_argument("-x", "--nx" , default=inputs['nx'] , type=int, help="Horizontal resolution")
    parser.add_argument( "-y", "--ny",  default=inputs['ny'],  type=int, help="Vertical resolution")
    parser.add_argument("--tmin", default=inputs['tmin'],  type=float, help=" Beginning of reconstruction")
    parser.add_argument("--tmax", default=inputs['tmax'], type=float, help="End of reconstruction" )
    parser.add_argument("-w", "--wrong_dets" ,  nargs='*', default=inputs['wrong_dets'],metavar='WRONG DET',
                        dest='wrong_dets', help="List of wrong detectors that are ignored, starts from 0 !")
    parser.add_argument("-d",  "--data_smooth", default=inputs['data_smooth'],type=float, dest='data_smooth'
                        ,help="Use moving average for input data over `n` bins" )
    parser.add_argument("-R", "--reconstruct" , default=inputs['reconstruct'],  dest='reconstruct',
                        type=int, help='Perform tomographic reconstruction (0/1)' )
    parser.add_argument("-P", "--preview",  default=inputs['show_preview'], dest='show_preview',
                        action='store_true', help="Show preview and quit")
    parser.add_argument("-G", "--gui", default=inputs['gui'], action='store_true' , 
                        help='Start graphical interface')
    parser.add_argument("-D", "--deriv", default=inputs['deriv'], action='store_true' , 
                        help='Plot graph of derivation')
    parser.add_argument("-r","--proj", default=inputs['plot_projection'], action='store_true' , 
                        help='Plot the projection space')
    
    parser.add_argument("-a", "--danis",default=inputs['danis'], type=float, 
                        help="Ratio of anisotropic matrices, B = sigmoid(-n)*Bperp+sigmoid(n)*Bpar "+\
                            " => The bigger number the more will reconstruction follow "+\
                                "the field For MFI is recommended ~ 4"  )
    parser.add_argument("-e", "--errorscale",default=inputs['error_scale'], dest='error_scale',
                        type=float , help=" Scale the expected errors, The higher number"+\
                            " the smoother reconstruction. Original value is 1")
    parser.add_argument("-p", "--sparallel", default=inputs['solid_parallel'], dest='solid_parallel'
                        ,action='store_true', help="Use faster parallel solve, but "+\
                            "unbreakable and without progressbar in GUI" )
    parser.add_argument("--postprocessing", default=inputs['postprocessing'], dest='postprocessing'
                        ,action='store_true', help='Perform postprocessing -- '+\
                            'total emissivity, position of center of emissivity, entropy, smoothness'  )
    parser.add_argument("--asymmetry", default=inputs['asymmetry'], dest='asymmetry',action='store_true',
                        help='Compute asymmetries of the SXR profile '  )
    parser.add_argument("--plot_svd",action='store_true' , default=inputs['plot_svd'] ,
                        help='Substract first mode of SVD from plotted results')
    parser.add_argument('--plot_poloidal_spectrum', action='store_true' , default=inputs['plot_poloidal_spectrum'],
                        help='Plot a poloidal mode spectrum as function of the radius')
    parser.add_argument('-u', '--undersampling',  dest="data_undersampling",
                        default=inputs['data_undersampling'], help=" Use every `n`-th snapshot" , type=int)
    parser.add_argument( '--sx',  default=0, 
                        help="Shift of magnetic field in X", type=float )
    parser.add_argument( '--sy',  default=0, 
                        help="Shift of magnetic field in Y" , type=float)
    parser.add_argument("--rem_back",  default=inputs['rem_back'], dest='rem_back',
                        action='store_true', help="Remove background from the reconstruction")
    parser.add_argument("--ifishmax", dest='ifishmax', default=inputs['ifishmax'],
                        type=int , help="Maximal number of MFI iterations")
    parser.add_argument( '--transform_order_a', default=inputs['transform_order_a'],type=int,
                        dest='transform_order_a', help="Angular order of transform basic vectors")
    parser.add_argument( '--transform_order_r', default=inputs['transform_order_r'],dest='transform_order_r',
                        help="Radial order of transform basic vectors" , type=int)
    parser.add_argument('-T', '--transform', default=inputs['transform_index'],dest='transform_index',
                        help="""Select orhogonal transformation - 0 - None, 1 - Abel, 2 - Cormack, 
                        3 - Fourier-Bessel""", type=int)
                        #, 4 - zoom, 5- artificial, 6 - in-out asymmetry, 7 - even polynom""" , type=int)
    parser.add_argument('--prune_dets', default=inputs['prune_dets'],dest='prune_dets',
                        help="Use only each n-th detector (for cameras)" , type=int)
    parser.add_argument('--plot_loglin', default=inputs['plot_loglin'],dest='plot_loglin',
                        action='store_true', help="Plot resultes in lin/log scale" )
    parser.add_argument('--sawtooths', default=inputs['sawtooths'],dest='sawtooths',
                        help="Use sawtooth detection algorithm" , action='store_true')
    parser.add_argument('--output_type', default=inputs['output_type'],dest='output_type',
                        help="file type of the plots when output is enable" , type=str)
    parser.add_argument('--positive_constrain', default=inputs['positive_constrain'],action='store_true',
                        dest='positive_constrain', help="force solution to be positive '+\
                        '(not compatible with solver 0 and 1)" )
    parser.add_argument('--impurities', default=inputs['impurities'],action='store_true',
                    dest='impurities', help="Impurity analysis from SXR" )
    
    parser.add_argument( "--dpi" , default=inputs['dpi'] , type=int, help="Resolution of the figures")
    parser.add_argument( "--img_size",  default=inputs['img_size'],  type=int, help="Size of the figures")
    parser.add_argument( "--tmp_folder",  default=inputs['tmp_folder'],  type=str, help="Temporal folder")

    args = parser.parse_args()

    for a in vars(args):
        inputs[a] = getattr(args,a)

    if inputs['error_scale'] <= 0: inputs['error_scale'] = 1e-2
 

    config.DEBUG = inputs['DEBUG']
    print('tmp,  ', args.tmp_folder)
    if args.DEBUG:
        rapid_blocks = 1

    wrong_dets =  inputs['wrong_dets']
    if wrong_dets != [] and isinstance(wrong_dets[0],str):
        wrong_dets = eval('r_['+''.join(wrong_dets)+']')-1
 
    config.wrong_dets_pref = []
    config.wrong_dets_defined = wrong_dets

    config.magfield_shift = args.sx,args.sy
  
    assert args.tmin <= args.tmax, "tmin is greater than tmax " + str(args.tmin) + '  ' + str(args.tmax)

    tok = loaddata(inputs)

    if not inputs['gui'] or inputs['reconstruct']:
        tok.prepare_tokamak()
        
        
    if hasattr(tok,'min_tvec'):
       inputs['tmin'] = max(tok.min_tvec,inputs['tmin'])
    if hasattr(tok,'max_tvec'):
       inputs['tmax'] = min(tok.max_tvec,inputs['tmax'])

    #print(tok.max_tvec)
    #print(inputs['tmin'], inputs['tmax'])
    #exit()

    if inputs['show_preview']:
        from matplotlib.backends.backend_qt5agg import FigureCanvasQT as FigureCanvas
        import matplotlib.pylab as plt
        preview(plt.gcf(),inputs, tok,True, True)
        plt.show()
        print('preview')
        plt.pause(100)

        sys.exit(0)

    elif inputs['deriv']:
        from graph_derivation import graph_derivation
        graph_derivation(4, tok)
        show()
        
        sys.exit(0)
        
    elif inputs['proj']:
        from geom_mat_setting import loadgeometry,plot_projection_space
        
        xchords, ychords, distance, nl,virt_chord  = loadgeometry(tok.geometry_path, list(tok.detectors_dict.keys()),100)   
        plot_projection_space(tok,xchords, ychords, virt_chord,nl)
        sys.exit(0)
        
    elif inputs['gui']:
        startGUI(inputs, tok)
        
    elif inputs['reconstruct'] or inputs['ntm']:
        startCUI(inputs, tok)

    print("""\n\n\nPlease, acknowledge this tomographic code by citing these papers: 
    
    Odstrcil, T., et al. "Optimized tomography methods for plasma 
        emissivity reconstruction at the ASDEX  Upgrade tokamak.
        " Review of Scientific Instruments 87.12 (2016): 123505.
    
    Odstrcil, M., et al. "Modern numerical methods for plasma 
        tomography optimisation."  Nuclear Instruments and Methods 
        in Physics Research Section A: Accelerators, Spectrometers,
        Detectors and Associated Equipment 686 (2012): 156-161.

    Thank you. 
    """)


def startGUI(inputs, tok):
    """   Start graphical interface: load data for choosen tokamak and open GUI.  In case of errors, try to load commandline interface

    :param inputs dict:  dictionary with all setting generated by `main` and `loadSetting`
    :raises ImportError: missing basic dependencies
    """



    if sys.platform == 'linux2':
        f=os.popen("xset -q")
        if len(f.readlines()) < 10:
            print('X server not availible\n Try to set export DISPLAY=":0.0" \n Switching to noGUI version')
            tok.prepare_tokamak()

            startCUI(inputs, tok)
            sys.exit(0)
            
            
            
            
    try:
        import PyQt5.QtCore,  PyQt5.QtGui ,PyQt5.QtWidgets
        import GUI
    except ImportError as details:
        try:
            import PyQt4.QtCore,PyQt4.QtGui 
            import GUI

        except ImportError as details:
            print('Install PyQt4 or PyQt5!!!!  Fallback to noGUI version')
            tok.prepare_tokamak()
            startCUI(inputs, tok)
            return
    try:
        GUI.main(inputs, tok)
    except:
        print("Loading GUI failed, try to remove cache (`cfg_file`.npz)")
        raise


def startCUI(inputs, tok):
    """
    Start command line user interface, load data for choosen tokamak and perform tomographic reconstruction.

    :param inputs dict:  dictionary with all setting generated by `main` and `loadSetting`
    """

    
    tmp_folder   = os.path.expanduser(os.path.expandvars(inputs['tmp_folder']))
    output_path  = os.path.expanduser(os.path.expandvars(inputs['output_path']))

    if not os.path.isdir( tmp_folder):
        try:
            os.mkdir(tmp_folder)
        except:
            print('error:os.mkdir(tmp_folder)')
            raise
    if not os.path.isdir( output_path):
        try:
            os.mkdir(output_path)
        except:
            print('error:os.mkdir(output_path)')
            raise
     
    
    inputs['postprocessing'] |= inputs['impurities']
    
    if inputs['reconstruct']:
        
        print('Reconstruction of shot '+str(inputs['shot'])+' started')
        output = tomography(inputs, tok)
        from make_graphs import make_graphs

        make_graphs(output)

        if inputs['plot_svd']:
            make_graphs(output,True)
                        
        if inputs['postprocessing'] and inputs['tsteps']>1:
            from make_graphs import postprocessing_plot
            postprocessing_plot(output)
            
        if inputs['plot_poloidal_spectrum'] and inputs['tsteps']>1:
            from make_graphs import CalcPoloidalModeSpectrum
            CalcPoloidalModeSpectrum(output)
                        
        if inputs['sawtooths'] and inputs['tsteps']>1:
            from sawtooths import sawtooths_detection
            sawtooths_detection(output)

        if inputs['impurities'] :
            from imp_analysis import imp_analysis
            imp_analysis(output)
        
        if inputs['plot_svd'] :
            from make_graphs import make_svd
            make_svd(*output)
        
        if inputs['asymmetry']:
            from asymmetries import  EvalAsymmetry
            EvalAsymmetry(  output)

            

            

if __name__ == "__main__":
    main()
