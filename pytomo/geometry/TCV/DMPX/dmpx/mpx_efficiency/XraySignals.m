function out=XraySignals(det)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function which calculates the efficiency (percentage of incoming energy detected) of a 
%  detector, as a function of the energy of the incoming photons
%  
%  SYNTAX
%
%  out=XraySignals(det)
%
%  INPUT
%
%  det	Detector string:
%	Everything before $ symbol is a filter, everything after
%	is an absorber. Material symbol (case is not important) 
%	is followed by thickness in MICRONS.
%	
%	LIST of materials for filters:
%	Be = beryllium
%	Al = aluminium
%	Si = silicium
%	He = helium
%	KrCH = 90% krypton + 10% methan
%	ArCH = 90% argon + 10% methan
%	XeCH = 90% xenon + 10% methan
%	AIR = air (mixture 78.1% N, 21% O, 0.015% C, 0.47% Ar)
%	
%	LIST of materials for absorbers:
%	
%	KrCHabs = 90% krypton + 10% methan
%	ArCHabs = 90% argon + 10% methan
%	XeCHabs = 90% xenon + 10% methan
%
%  EXAMPLES
%  
%  * two filters system with 5000 micron Ar detector  
%  det = 'al200, Be200 $ Ar5000'
%	 
%  * DMPX top detector with KrCH absorber
%  det = 'He229600, Be100 $ KrCHabs8000'
%	 
%  * DMPX bottom detector with KrCH absorber
%  det = 'He229600, Be200, KrCH8000, AIR10000 $ KrCHabs7600'
%
%
%  OUTPUT
%   
%  out.response=out.tr.*out.abso	system efficiency 
%  out.abso				absorber efficency
%  out.tr				filter transmision 
%  out.ev				energy vector used
% 
%  From Alexei Zabolotsky routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    out=[];
    
    if nargin<1 
        help XraySignals.m
        return
    end    
 
	[filt,abso,error]=XraySignals_Strings(det);
	
	if error==1
       disp('XraySignals ERROR: cannot process detectors string') 
       return
	end   
	
    e=linspace(200,2e5,3000);   % energy vector (eV)

%***************  get filter transmission characteristics *********************
    for j=1:length(filt)
         [tx,error]=XraySignals_SxCurves(filt(j).material,filt(j).thick,e);
         tr(j,:)=tx;      % each filter transmission as a function of frequency
	end
	tr=prod(tr,1);  % total filter transmission as a function of frequency

%***************  get detector absorbtion characteristics *********************
	
    for j=1:length(abso)
         [tx,error]=XraySignals_SxCurves(abso(j).material,abso(j).thick,e);
         diode_eff(j,:)=(1-tx);  % each diode efficency as a function of frequency
	end
	diode_eff=prod(diode_eff,1); % total diode efficency as a function of frequency


%****************************** Form output ***************************************    

    out.abso=diode_eff;     	 % diode efficency
    out.tr=tr;                    % filter transmision 
    out.ev=e;                     % energy vector used
    out.response=tr.*diode_eff;   % system energy respond 
