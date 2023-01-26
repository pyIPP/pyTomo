

## Install:
    module load git
    git clone https://github.com/pyIPP/tomo.git
    cd ./pytomo
    ./pytomo.sh



## Description:
pyTomo is advanced code for tomographical inversion from various line integrated diagnostics like SXR or bolometers on DIII-D. The main advantages are a simple GUI, very fast inversion and a high accuracy (depends on the quality of the input data).

 
## GUI
Run GUI by module load pytomo
and pytomo
If you want to open it directly with some discharge loaded, use

    pytomo -G shot

### Load data
First, select in Input panel discharge, and source and press Data Setting. In the opened window can be selected corrupted channels by left double-click (place back by right double-click) or by specifying their numbers. Also, you should select the time range, downsampling, and smoothing, too many timeslices can make server administrators angry!

In Output panel can be specified where to store the results (if any), set publication quality plots (slow!!), show separate plot for each timeslice (can be also slow) etc... Other panels specify details of the tomography inversion, default options are usually fine for DIII-D SXR. For details, check the presentation lower on this page.

The last step is pressing of START button. In case of any problems, contact the code developer.

 ![1200px-Pytomo](https://user-images.githubusercontent.com/32073690/214947537-5ddfe8a3-85c8-4e33-adaa-c1573b7a925a.png)

### Command line
Everything what can be set in GUI can be defined also by command line switches, write:

    /fusion/projects/codes/pytomo/pytomo.sh -h
    
for a list of availible options. Here is an example for DIII-D bolometers

    /fusion/projects/codes/pytomo/pytomo.sh 163303 -t 23 --lambda_solver 0.2 -u 5 -d 5 --noboundary -S 7 --positive_constrain -G -B 100

### PowerPoint manual
Manual for AUG, but most of it is valid for DIII-D as well. File:SXR tomography on AUG.pptx


## References
 Odstrcil, T., et al. "Optimized tomography methods for plasma emissivity reconstruction at the ASDEX Upgrade tokamak. " Review of Scientific Instruments 87.12 (2016): 123505.
 
 Odstrcil, M., et al. "Modern numerical methods for plasma tomography optimisation." Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment 686 (2012): 156-161.
