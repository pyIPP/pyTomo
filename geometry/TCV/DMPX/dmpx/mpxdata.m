% MPXDATA allows you to load the DMPX calibrated signal, detectors high voltage, geometry,
% filters, an estimation of rhopsi and rhovol corresponding to the lines-of-sight as well as
% the detectors efficiency. MPXDATA works only for shots > 20030.
%
% SYNTAX
% [mpx] = mpxdata(shot,action,varargin)
%
% INPUTS
% shot		= Shot number
% action    = String containing the key character for the actions you want to execute. The actions are:
%       's':	Loads the MPX calibrated signal.
%       Channels are sorted from HFS to LFS.
%       Warning:
%       * For shots < 27127, the calibration does not include the voltage applied to the detector.
%       * For shots > 27128, it does and the signal of two shots with different voltage can be compared.
%       'v':	Loads the detector high voltage value [V] and corresponding gain.
%       'g':	Loads the detector geometry (slit and wires position from HFS to LFS) [m].
%               Creates also a dcd object to be used with the psitbx.
%       'f':	Loads the filters in front of the detectors (thicknesses in [m]).
%               + the position of the mobile absorber holder
%               + the gas in the detector
%               Use the 'pos' option to force the position of the mobile absorber holder.
%       'r':	Calculates the rho_psi to which the diagnostic chords are tangent (psitbx times), from HFS to LFS.
%               Use the 'rtime', 'liuqe', 'psi', 'opt' and 'vol' options to specify the rho parameters.
%       'e':	Loads the detector efficiency (probability that the detection gas absorbs an incoming photon)
%               as a function of the energy of incoming photons. The MPX output signal is proportional to
%               the detector efficiency if the escape peak is neglected.
% varargin	= Parameters for the choosen actions. You can specify:
%       'time'		Time interval (ex: [0.1 0.9]) for the signal. Default = full shot.
%		'freq'		Acquisition frequency: 
%				0: slow 20kHz
%				1: fast 200kHz (if available, else slow)
%				2: ask if both are available, else slow (default)
%		'detec'		Detector selection:
%				'top':	top detector
%				'bot':	bottom detector
%				'both':	both detectors (default)
%		'chords'	list of chords to consider (ex: [32 48 63]). Default = all chords considered.
%				[1:64] -> top detector
%				[65:96] -> bottom detector
%				NB: The 'chords' option overrules the 'detec' option.
%		'pos'		To force the position of the mobile absorber holder for the efficiency calculation.
%		'rtime'		Time interval for the rho_psi calculation (ex: [0.1 0.9]). It then calculates
%				a value of rho on each psitbx time. Default = value used for action 's'.
%		'liuqe'		Liuqe version (1, 2, 3) for the rho_psi calculatin. Default = 1.
%		'psi'		Flux description in cylindrical coordinates ('01' psitbxpsi object) for the rho_psi
%				calculation. Optional, useful if you already have psi loaded in your workspace.
%				Ex: psi=psitbxtcv(shot,time,'01').
%		'opt'		Option for the rho_psi calculation
%				-1: to get rho from -1 (HFS) to 1 (LFS). Default.
%				+1: to get rho from 0 to 1.
%		'vol'		Non zero value to get also rho_vol. Default = 0.
%
% OUTPUT
% mpx		= structure containing the desired informations
%
% EXAMPLES
% [mpx] = mpxdata(27991,'svgf')
% [mpx] = mpxdata(27991,'sr','time',[0.9 1.1],'freq',1,'detec','bot','vol',1)
% [mpx] = mpxdata(24848,'g','detec','top')
% [mpx] = mpxdata(36989,'s','time',[0.9 0.91],'chords',[1 32 36 65 94])
%
% Written by Y. Camenen, CRPP, EPFL.
% Modified by L. Curchod, CRPP, EPFL:	15.07.2008	Correction of the error when loading rho_psi between -1 and 1.
% Modified by L. Curchod, CRPP, EPFL:	14.08.2008	Use a new mpx_calib_dec05_bis.mat file calculated with data for 2450 V as well. 
% Modified by L. Curchod, CRPP, EPFL:	15.08.2008	Corrected to sort the gain coefficients for the oct. 2005 calibration case.
% Modified by L. Curchod, CRPP, EPFL:	03.03.2009	New path to the calibration files on PC70
% Modified by Y. Camenen, Warwick, UK:	10.07.2009	Add option to load selected chords
% Modified by J. Kamleitner, CRPP, EPFL:    3.2013  moved calibration files to crpptbx to be independent of user home directories

function [mpx]=mpxdata(shot,action,varargin)


%
% Check for non-stored shots
%
if (exist('shot','var') && ~isnumeric(shot))
	error('First argument must be the shot number.')
elseif shot<20030
	error('This program works only for shot >= #20030.')
elseif (shot>=34988)&(shot<=35438)
	warning('DMPX signal has been lost for shots #34988 to #35428.');
	action(strfind(action,'s'))=''; % Withdraw the action of signal loading;
elseif (shot>=39113)&(shot<=39123)
	%warning('DMPX signal has been lost for shots #39113 to #39123.');
	%action(strfind(action,'s'))=''; % Withdraw the action of signal loading;
end

fprintf('\n+ MPXDATA for shot %d.\n',shot);

%
% Define defaults values
%
mpxpath='/home/matlab/crpptbx-7.6.0/dmpx/';
def.action = 's';
def.time='**';
def.rtime=[0 10];
def.det='both';
def.freq=2; % (0=slow, 1=fast if available, 2=ask, if fast is available)
def.chords=[1:96]; 
def.liuqe=1;
def.opt=-1;
def.vol=0;
def.top.nchords=64;
def.bot.nchords=32;
abso.material={'Al','Be','Be','nothing'};
abso.thickness=[308 550 125 NaN]*1e-6;
%
% Fill in the variables with the values in varargin
%
if nargin > 2
    if rem(nargin,2)==0,
        for ii=1:2:nargin-2,
            eval([varargin{ii} '=varargin{ii+1};'])
        end
    else
        error('mpxdata:WrongInput','Wrong number of input arguments')
    end
end
%
% Default values
%
if ~exist('detec','var')
	detec = def.det;
end
if ~exist('action','var')
    action = def.action;
end
%
% Check 'chord' and 'detec' consistency
%
if exist('chords','var') % If 'chords' is specified, it overrules 'detec'
	if exist('detec','var')
		old_detec = detec;
	end
	if any(chords<1)||any(chords>def.top.nchords+def.bot.nchords)
		error(['Option ''chords'' must only include numbers between 1 and ' num2str(def.top.nchords+def.bot.nchords) ])
	end
	istop=any(chords<=def.top.nchords);
	isbot=any(chords>def.top.nchords);
	switch istop+isbot*2
	case 1
		detec='top';
	case 2
		detec='bot';
	case 3    	
		detec='both';
	end
	if (exist('old_detec','var') && ~strcmp(old_detec,detec))
		warning(['Option ''detec'' has been overruled from ''',old_detec,''' to ''',detec,''' by the ''chords'' selection.']);
	end
else
    if strcmp(detec,'both')
        chords=def.chords;
    elseif (strcmp(detec,'top'))
        chords=1:def.top.nchords;
    elseif (strcmp(detec,'bot'))
        chords=(1:def.bot.nchords)+def.top.nchords;
    end
end
if shot<26609
	if strcmp(detec,'both'), 
		if shot<26555
			fprintf('- Load data for top detector only: Duplex MPX not installed for shot < #26555.\n');
		else
			fprintf('- Load data for top detector only: No bottom detector HV source for shot < #26609.\n');
		end
		detec='top';
	elseif strcmp(detec,'bot')
		if shot<26555
			error('No bottom detector: Duplex MPX not installed for shot < #26555.')
		else
			error('No bottom detector: No bottom detector HV source for shot < #26609.')
		end
	end
end
%
% To load the MPX geometry before computing rho
%
i_r=strfind(action,'r');
if ~isempty(i_r),
	i_g=strfind(action,'g');
	if isempty(i_g),
		action(end+1)='g';
		i_g=length(action);
	end
	if i_g>i_r
		action(i_r)='g';
		action(i_g)='r';
	end
end
%
% To load the MPX high voltage value before calibrating the signal
%
i_s=strfind(action,'s');
if ~isempty(i_s),
	i_v=strfind(action,'v');
	if isempty(i_v),
		action(end+1)='v';
		i_v=length(action);
	end
	if i_v>i_s
		action(i_s)='v';
		action(i_v)='s';
	end
end
%
% To load the MPX filters before computing the energy response
%
i_e=strfind(action,'e');
if ~isempty(i_e),
	i_f=strfind(action,'f');
	if isempty(i_f),
		action(end+1)='f';
		i_f=length(action);
	end
	if i_f>i_e
		action(i_e)='f';
		action(i_f)='e';
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     LOAD DATA                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:length(action)
switch action(ii)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Load the DMPX calibrated signal       %	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
case 's' 
    fprintf('\ncase s\n');
mdsopen(shot);
mdsdata('reset_public()');
%	 
% Some initalisations and tests
%
fprintf('- Load the MPX calibrated signal.\n');
if exist('freq','var')~=1
	freq=def.freq;
end
if exist('time','var')==1, 
	if shot<20939 % time base was incorrect
		time(1)=(time(1)+0.04)/0.9697-0.04;
		time(2)=(time(2)+0.04)/0.9697-0.04;
	end
else
	time=def.time;
end
echant=['[' num2str(time(1)) ':' num2str(time(2)) ':*]'];
if( shot<24087 || shot>24725 )
	DTNE1='\ATLAS::DT100_NORTHEAST_001:CHANNEL_0';
	DTNE2='\ATLAS::DT100_NORTHEAST_002:CHANNEL_0';
	DTNE3='\ATLAS::DT100_NORTHEAST_003:CHANNEL_0';
	DTNE1_fast='\ATLAS::DT100_NORTHEAST_001:SELECTED:CHANNEL_0';
	DTNE2_fast='\ATLAS::DT100_NORTHEAST_002:SELECTED:CHANNEL_0';
	DTNE3_fast='\ATLAS::DT100_NORTHEAST_003:SELECTED:CHANNEL_0';
else
	DTNE1='\ATLAS::DT100_SOUTHWEST_001:CHANNEL_0';
	DTNE2='\ATLAS::DT100_SOUTHWEST_002:CHANNEL_0';
	DTNE1_fast='\ATLAS::DT100_SOUTHWEST_001:SELECTED:CHANNEL_0';
	DTNE2_fast='\ATLAS::DT100_SOUTHWEST_002:SELECTED:CHANNEL_0';
end
%     
% Test the DTACQ trigering mode
%
mode1=mdsdata([DTNE1(1:end-9) 'MODE']);
mode2=mdsdata([DTNE2(1:end-9) 'MODE']);
if(strcmp(detec,'bot') || strcmp(detec,'both')) 
	mode3=mdsdata([DTNE3(1:end-9) 'MODE']); 
else
	mode3=4;
end;
mode=mode1*mode2*mode3; %mode=64 ->ok
% if shot>=24087&shot<25095
if shot>=24087&mode~=64
	warning('DTACQ not in mode 4')
	disp('Random temporal gap (5 to 25 us) between the two or three MPX acquisition cards.')
	disp('Random temporal gap (0.1 to 0.2ms) between DTACQ and TCV.')
end
%
% Acquisition frequency
%
if (shot<34988) % After TCV BIG OPENING 2007-2008, fast data are in standard nodes
	if freq~=0
		test_fast=mdsdata(['node_size("\' DTNE1_fast '01")']); %test whether the first channel of the first card has fast acquisition
        if(isnumeric(test_fast))
            if(test_fast>0)
                flag1=0;
                if freq==2,
                    freq=input(['Do you want standard acquisition f = 20kHz (0) or fast acquisition at f = 200 kHz (1) ?']);
                end
            else
                flag1=1;
            end
        else
            flag1=1;
        end
        if(flag1)
			disp('Only standard acquisition (f=20kHz) available')
			freq=0;
        end
	end
else % After TCV BIG OPENING 2007-2008, fast data is always available in standard nodes
	if freq==2,
		freq=input(['Do you want decimated data at f = 20kHz (0) or full data at f = 200 kHz (1) ?']);
	end
end
%	 
% Load, sort and calibrate TOP data
%
if strcmp(detec,'top')|strcmp(detec,'both'),
	if freq==1,
		if (shot<34988) % Before TCV BIG OPENING 2007-2008, fast data are in selected nodes
			DTNE1=DTNE1_fast;
			DTNE2=DTNE2_fast;
		end
	else
		if (shot>=34988) % After TCV BIG OPENING 2007-2008, fast data are in standard nodes
			echant=['[' num2str(time(1)) ':' num2str(time(2)) ':0.00005]']; % Need to decimate for 20 kHz
		end
	end 
	%
	% Load timebase and check if node is full
	%
    datatest = tdi([DTNE1,'01',echant]);
	if isempty(datatest.data)
		warning('No TOP detector signal.');
		load_opt = 0;
	else
		mpx.top.signal.dim{1} = datatest.dim{1};
		mpx.top.signal.dim{2} = [1:def.top.nchords];
		load_opt = 1;
	end
	if load_opt
		if (shot>19923 & shot<20939)  % time base incorrect, need to rescale
			mpx.top.signal.dim{1}=(mpx.top.signal.dim{1}+0.04)*0.9697-0.04;
		end;
		l_t=length(mpx.top.signal.dim{1});
		tmpdata=NaN.*zeros(l_t,length(mpx.top.signal.dim{2}));
		II=1:32;
		ch=int2str(II');
		ch(ch==' ')='0';	
		%
		% Sort the channels from HFS to LFS
		%
		if shot<26555 %MPX, calib_coeff is already sorted for these shots
			II(1:2:63)=[33:64];
			II(2:2:64)=[32:-1:1];
		elseif shot>=26555 %DMPX, calib_coeff is not sorted for these shots
			II(1:2:63)=[33:64];
			II(2:2:64)=[16:-1:1 32:-1:17];
			calib_coeff_t=calib_coeff_t(II);	% Ordering calib_coeff for the top detector
			mpx.top.gain.C=mpx.top.gain.C(II);	% Ordering gains for the top detector (useful only for calibration of oct. 05
		end
		%
		% Load the raw data, not sorted, not calibrated   
		% NB: only load the data for the channels specified in 'chords'
		%
        
        for iii=1:32,
            if any(find(II-iii==0)-chords==0)
                %[DTNE1 ch(ii,:) echant]
                tmpdata(:,iii)=mdsdata([DTNE1 ch(iii,:) echant]);
            end
            if any(find(II-iii-32==0)-chords==0)
                tmpdata(:,iii+32)=mdsdata([DTNE2 ch(iii,:) echant]); 
            end	
        end 
        
        %
		% Calibrate
		%
		mpx.top.signal.data=tmpdata(:,II).*repmat(calib_coeff_t,l_t,1)./repmat(mpx.top.gain.C,l_t,1);
		%
		% Some special cases
		%
		if shot>=25102 && shot<25933, 
			mpx.top.signal.data(:,36)=NaN;
			disp('Channel 36 was missing for this shot') 
		end	
		if shot>=27127 && shot<=28124 
			mpx.top.signal.data(:,[4 64 62 60 58 56])=NaN; %missing channels, problem with cable 2, repaired by DF in Dec 2004
			if shot>=27185 %one more channel missing !...
				mpx.top.signal.data(:,[44])=NaN;
			end 
		end 
		if shot>=28219 && shot<31446, 
			mpx.top.signal.data(:,[19 21])=NaN;
			disp('Channel 19 and 21 were missing for this shot') 
		end
	end
end	
%	 
% Load and calibrate BOTTOM data (no need to sort)
%
if strcmp(detec,'bot') || strcmp(detec,'both'),
	if freq==1,
		if (shot<34988) % Before TCV BIG OPENING 2007-2008, fast data are in selected nodes
			DTNE3=DTNE3_fast;
		end
	else
		if (shot>=34988) % After TCV BIG OPENING 2007-2008, fast data are in standard nodes
			echant=['[' num2str(time(1)) ':' num2str(time(2)) ':0.00005]']; % Need to decimate for 20 kHz
		end
	end 
	datatest = tdi([DTNE3,'01',echant]);
	if isempty(datatest.data)
		warning('No BOTTOM detector signal.');
		load_opt = 0;
	else
		mpx.bot.signal.dim{1} = datatest.dim{1};
		mpx.bot.signal.dim{2} = [1:def.bot.nchords];
		load_opt = 1;
	end
	if load_opt
		l_t=length(mpx.bot.signal.dim{1});
		tmpdata=NaN.*zeros(l_t,length(mpx.bot.signal.dim{2}));
		II=1:32;
		ch=int2str(II');
		ch(ch==' ')='0';	
		%
		% Load the raw data, sorted from HFS to LFS, non calibrated   
		%
		for iii=1:32,
			if any(def.top.nchords+find(II-iii==0)-chords==0)
				tmpdata(:,iii)=mdsdata([DTNE3 ch(iii,:) echant]);
			end
		end
		% 
		% Calibrate
		%
		mpx.bot.signal.data=tmpdata(:,II).*repmat(calib_coeff_b,l_t,1)./repmat(mpx.bot.gain.C,l_t,1);
		%
		% Some special cases
		%
		if shot>=30759 && shot<31446, 
			mpx.bot.signal.data(:,24)=NaN;
			disp('Calibration problem for channel 24.') 
		end
	end
end
%	 
% merge top and bottom data
% the bottom data is therefore rescaled by R(2kV)/2=0.2848 according to
% coeff.m
%
if(strcmp(detec,'top') || strcmp(detec,'both'))
    mpx.signal.dim{1}=mpx.top.signal.dim{1};
else
    mpx.signal.dim{1}=mpx.bot.signal.dim{1};
end
mpx.signal.dim{2}=chords;
if(strcmp(detec,'both'))
    mpx.signal.data=[mpx.top.signal.data(:,chords(chords<=def.top.nchords)) mpx.bot.signal.data(:,chords(chords>def.top.nchords)-def.top.nchords)];
elseif(strcmp(detec,'top'))
    mpx.signal.data=mpx.top.signal.data(:,chords(chords<=def.top.nchords));
else
    mpx.signal.data=mpx.bot.signal.data(:,chords(chords>def.top.nchords)-def.top.nchords);
end

mdsclose;
% fprintf('- MDS disconnect.\n');
% mdsdisconnect;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Load the detectors high voltage value       %	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'v' 	
    fprintf('\ncase v\n');
mdsopen(shot);
mdsdata('reset_public()');
fprintf('- Load the detectors high voltage value\n');
if shot<26765 %After #26575, the real voltage is indicated in the Vista window, before it was the reference voltage
	mm=500;
else
	mm=1;
end 
if strcmp(detec,'top') || strcmp(detec,'both'),
	mpx.top.voltage=mdsdata('\VSYSTEM::TCV_PUBLICDB_R["ILOT:DAC_V_01"]')*mm;
 	if shot==32035 %special case were mpx voltage was switched of before the acquisition was finished
		mpx.top.voltage=2000;
	end
	mpx.top.gain.C=zeros(1,def.top.nchords);
	I=find(shot>=[20030 23323 26555 27127 29921 30759 31446]); % Better than "if...else...end": the array can grow
	switch I(end)
	case 1
		disp('Detector gain dependence on the high voltage value not included in the signal calibration');
		load(sprintf('%smpx_calib_first.mat',[mpxpath 'calibration_used/']),'calib_coeff');
		calib_coeff_t=calib_coeff;
		mpx.top.gain.C(:)=1;
	case 2
		disp('Detector gain dependence on the high voltage value not included in the signal calibration');
		load(sprintf('%smpx_calib_sept02.mat',[mpxpath 'calibration_used/']),'calib_coeff');
		calib_coeff_t=calib_coeff;
		mpx.top.gain.C(:)=1;
	case 3
        	disp('Detector gain dependence on the high voltage value not included in the signal calibration');
        	warning('There were leaks in the top detector wire chamber for 26554<shot<27128')
		disp('Calibration is not very meaningful')
		load(sprintf('%smpx_calib_may04.mat',[mpxpath 'calibration_used/']),'calib_coeff')
		calib_coeff_t=calib_coeff;
		mpx.top.gain.C(:)=1;	       
	case 4
        	disp('Same gain dependence on the high voltage value taken for each channel')
		
		load(sprintf('%smpx_calib_july04.mat',[mpxpath 'calibration_used/']),'calib_coeff','R')
		if strcmp(getenv('USER'),'goodman')|strcmp(getenv('USER'),'camenen')|strcmp(getenv('USER'),'curchod'),
			ch=1;
			ch=input('Hello Tim, do you want to make NaN disapear in the MPX profiles (1/0) ?');
			if ch==1
				c_old=calib_coeff;
				load(sprintf('%smpx_calib_may05.mat',[mpxpath 'calibration_used/']),'calib_coeff')      
				c_new=calib_coeff;
				tmp=isnan(c_old);
				calib_coeff(tmp)=c_new(tmp);
				calib_coeff(~tmp)=c_old(~tmp); 
			end 
		end
		calib_coeff_t=calib_coeff;
		load(sprintf('%smpx_calib_may05.mat',[mpxpath 'calibration_used/']),'C','V')
		C=mean(C(:,1:def.top.nchords),2); % Use the same gain for each channel
		%tmp=interpos(13,V,log(C),[mpx.top.voltage mpx.top.voltage]);
		%mpx.top.gain.C(:)=exp(tmp(1));
        mpx.top.gain.C(:)=exp(spline(V,log(C),mpx.top.voltage));
		mpx.top.gain.R=R;
	case 5
		disp('Same gain dependence on the high voltage value taken for each channel')
		load(sprintf('%smpx_calib_may05.mat',[mpxpath 'calibration_used/']),'calib_coeff','C','V')
		calib_coeff_t=calib_coeff;
		C=mean(C(:,1:def.top.nchords),2); % Use the same gain for each channel
		%tmp=interpos(13,V,log(C),[mpx.top.voltage mpx.top.voltage]);
		%mpx.top.gain.C(:)=exp(tmp(1));
        mpx.top.gain.C(:)=exp(spline(V,log(C),mpx.top.voltage));
		load(sprintf('%smpx_calib_july04.mat',[mpxpath 'calibration_used/']),'R') %use the previous relative calibration
		mpx.top.gain.R=R;
	case 6
		%
		% In this case, the different behaviour of the wires is contained in the matrix of gains.
		% The calibration coefficients are in a vector: one value per wire, same value for all tensions.
		%
		disp('Gain dependence on the high voltage value calibrated for each channel')
		disp('Leaks in the bottom detector, no relative calibration of the two detectors')
		load(sprintf('%smpx_calib_oct05.mat',[mpxpath 'calibration_used/']),'calib_coeff','C','V')
		calib_coeff_t=calib_coeff;
		for jj=1:def.top.nchords % Interpolation to get the proper gains wrt to the high tension value
			%tmp=interpos(13,V,log(C(:,jj)),[mpx.top.voltage mpx.top.voltage]);
			%mpx.top.gain.C(jj)=exp(tmp(1));
            mpx.top.gain.C(jj)=spline(V,log(C(:,jj)),mpx.top.voltage);
		end
		mpx.top.gain.R=NaN;
	case 7
		%
		% In this case, the different behaviour of the wires is contained in the matrix of calibration coefficients.
		% The gains are in a vector: one value per tension, same value for all wires.
		%
		disp('Gain dependence on the high voltage value calibrated for each channel')
		load(sprintf('%smpx_calib_dec05_bis.mat',[mpxpath 'calibration_used/']),'calib_coeff_top','C_top_av','V_top','R');
		for jj=1:def.top.nchords, % Interpolation to get the proper calibration coefficient wrt the high tension value
			%tmp=interpos(13,V_top,calib_coeff_top(:,jj),[mpx.top.voltage mpx.top.voltage]);
			%calib_coeff_t(jj)=tmp(1);
            calib_coeff_t(jj)=spline(V_top,calib_coeff_top(:,jj),mpx.top.voltage);
		end
		%tmp=interpos(13,V_top,log(C_top_av),[mpx.top.voltage mpx.top.voltage],1e6);  % Interpolation to get the proper gains wrt to the high tension value
		%mpx.top.gain.C(:)=exp(tmp(1));
        mpx.top.gain.C(:)=exp(spline(V_top,log(C_top_av),mpx.top.voltage));
		mpx.top.gain.R=R;		
	end
end
if strcmp(detec,'bot') || strcmp(detec,'both'),
	mpx.bot.voltage=mdsdata('\VSYSTEM::TCV_PUBLICDB_R["ILOT:DAC_V_02"]')*mm;
	if shot==32035 %special case were mpx voltage was switched of before the acquisition was finished
		mpx.top.voltage=2000;
	end  
	mpx.bot.gain.C=zeros(1,def.bot.nchords);
	I=find(shot>=[26555 27127 29921 30759 31446]);  	
	switch I(end)
  	case 1
		disp('Detector gain dependence on the high voltage value not included in the signal calibration')
		warning('There were leaks in the bottom detector wire chamber for 26554<shot<27128')
		disp('Calibration is not very meaningful')
		load(sprintf('%smpx_calib_may04.mat',[mpxpath 'calibration_used/']),'calib_coeff')
		calib_coeff_b=calib_coeff(def.top.nchords+1:def.top.nchords+def.bot.nchords);
		mpx.bot.gain.C(:)=1;
	case 2
		disp('Same gain dependence on the high voltage value taken for each channel')
		load(sprintf('%smpx_calib_july04.mat',[mpxpath 'calibration_used/']),'calib_coeff','R') 
		calib_coeff_b=calib_coeff(def.top.nchords+1:def.top.nchords+def.bot.nchords);
		load(sprintf('%smpx_calib_may05.mat',[mpxpath 'calibration_used/']),'C','V')
		C=mean(C(:,def.top.nchords+1:def.top.nchords+def.bot.nchords),2); %use the same gain for each channel
		%tmp=interpos(13,V,log(C),[mpx.bot.voltage mpx.bot.voltage]);
		%mpx.bot.gain.C(:)=exp(tmp(1));
        mpx.bot.gain.C(:)=exp(spline(V,log(C),mpx.bot.voltage));
		mpx.bot.gain.R=R;
	case 3
		disp('Same gain dependence on the high voltage value taken for each channel')
		load(sprintf('%smpx_calib_may05.mat',[mpxpath 'calibration_used/']),'calib_coeff','C','V')
		calib_coeff_b=calib_coeff(def.top.nchords+1:def.top.nchords+def.bot.nchords);
		C=mean(C(:,def.top.nchords+1:def.top.nchords+def.bot.nchords),2); %use the same gain for each channel
		%tmp=interpos(13,V,log(C),[mpx.bot.voltage mpx.bot.voltage]);
		%mpx.bot.gain.C(:)=exp(tmp(1));
        mpx.bot.gain.C(:)=exp(spline(V,log(C),mpx.bot.voltage));
		load(sprintf('%smpx_calib_july04.mat',[mpxpath 'calibration_used/']),'R') %use the previous relative calibration
		mpx.bot.gain.R=R;
	case 4
		%
		% In this case, the different behaviour of the wires is contained in the matrix of gains.
		% The calibration coefficients are in a vector: one value per wire, same value for all tensions.
		%
		disp('Gain dependence on the high voltage value calibrated for each channel')
		disp('Leaks in the bottom detector, no relative calibration of the two detectors')
		load(sprintf('%smpx_calib_oct05.mat',[mpxpath 'calibration_used/']),'calib_coeff','C','V')
		calib_coeff_b=calib_coeff(def.top.nchords+1:def.top.nchords+def.bot.nchords);
		for jj=1:def.bot.nchords,
			%tmp=interpos(13,V,log(C(:,jj+def.top.nchords)),[mpx.bot.voltage mpx.bot.voltage]);
			%mpx.bot.gain.C(jj)=exp(tmp(1));
            mpx.bot.gain.C(jj)=exp(spline(V,log(C(:,jj+def.top.nchords)),mpx.bot.voltage));
		end
		mpx.bot.gain.R=NaN;
	case 5
		%
		% In this case, the different behaviour of the wires is contained in the matrix of calibration coefficients.
		% The gains are in a vector: one value per tension, same value for all wires.
		%
		disp('Gain dependence on the high voltage value calibrated for each channel')
		load(sprintf('%smpx_calib_dec05_bis.mat',[mpxpath 'calibration_used/']),'calib_coeff_bot','C_bot_av','V_bot','R')
		for jj=1:def.bot.nchords,
			%tmp=interpos(13,V_bot,calib_coeff_bot(:,jj),[mpx.bot.voltage mpx.bot.voltage]);
			%calib_coeff_b(jj)=tmp(1);
            calib_coeff_b(jj)=spline(V_bot,calib_coeff_bot(:,jj),mpx.bot.voltage);
		end
		%tmp=interpos(13,V_bot,log(C_bot_av),[mpx.bot.voltage mpx.bot.voltage],1e6);
		%mpx.bot.gain.C(:)=exp(tmp(1));    
        mpx.bot.gain.C(:)=exp(spline(V_bot,log(C_bot_av),mpx.bot.voltage));
		mpx.bot.gain.R=NaN;
	end
end
mdsclose;
% fprintf('- MDS disconnect.\n');
% mdsdisconnect;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Load the detectors geometry [m]       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'g'
    fprintf('\ncase g\n');
fprintf('- Load the detectors geometry.\n');
I=find(shot>=[20030 23806 26555]); % Better than "if...else...end": the array can grow
switch I(end) 
case 1
	disp('MPX geometry from shot 20030 to 23806')
	mpx.top.geom.wires.R=[0.943:-0.002:0.817];
	mpx.top.geom.wires.Z=repmat(-1.558,1,64);
	mpx.top.geom.wires.length=32e-3; %check it !
	mpx.top.geom.slit.R=0.88;
	mpx.top.geom.slit.Z=-1.1; 
	mpx.top.geom.slit.dR=2e-3; %slit width
	mpx.top.geom.slit.length=40e-3; %slit length (toroidal extension)
case 2
	disp('MPX geometry from shot 23807 to 26554')
	mpx.top.geom.wires.R=[0.943:-0.002:0.817];
	mpx.top.geom.wires.Z=repmat(-1.14921,1,64);
	mpx.top.geom.wires.length=32e-3; %check it !
	mpx.top.geom.slit.R=0.88;
	mpx.top.geom.slit.Z=-0.92261;
	mpx.top.geom.slit.dR=2e-3; %slit width
	mpx.top.geom.slit.length=32e-3; %slit length (toroidal extension)
case 3
	disp('DMPX geometry since shot 26555')
	mpx.top.geom.wires.R=[0.943:-0.002:0.817];
	mpx.top.geom.wires.Z=repmat(-1.15621,1,64);
	mpx.top.geom.wires.length=32e-3; %check it !
	mpx.top.geom.slit.R=0.88;
	mpx.top.geom.slit.Z=-0.92261;
	mpx.top.geom.slit.dR=2e-3; %slit width
	mpx.top.geom.slit.length=32e-3; %slit length (toroidal extension)
	mpx.bot.geom.wires.R=[0.942:-0.004:0.816];
	mpx.bot.geom.wires.Z=repmat(-1.17421,1,32);
	mpx.bot.geom.wires.length=32e-3; %check it !
	mpx.bot.geom.slit.R=0.88;
	mpx.bot.geom.slit.Z=-0.92261;
	mpx.bot.geom.slit.dR=2e-3; %slit width
	mpx.bot.geom.slit.length=32e-3; %slit length (toroidal extension)
end
%
% Create the dcd object
%
nd=1000;
if (strcmp(detec,'top') || strcmp(detec,'both'))
	rs=repmat(mpx.top.geom.slit.R,1,def.top.nchords); %slit position R
	zs=repmat(mpx.top.geom.slit.Z,1,def.top.nchords); %slit position Z
	rd=mpx.top.geom.wires.R; %wire position R
	zd=mpx.top.geom.wires.Z; %wire position Z
	phid=zeros(1,def.top.nchords); %toroidal location angle of wires
	tvd=zeros(1,def.top.nchords); %toroidal angles of lines of sight
	pvd=atan((zs-zd)./(rs-rd)); %poloidal angles of lines of sight
	pvd(pvd<0)=pi+pvd(pvd<0);
	pvd=pi-pvd; %angle between the horizontal plane and the chord, sens trigonometrique
	mpx.top.geom.dcd=psitbxdcd(rd,zd,phid,pvd,tvd,nd);
    if (strcmp(detec,'top'))
        rd=rd(chords);
        zd=zd(chords);
        phid=phid(chords);
        pvd=pvd(chords);
        tvd=tvd(chords);
    end
end
if (strcmp(detec,'bot') || strcmp(detec,'both'))
	rs=repmat(mpx.bot.geom.slit.R,1,def.bot.nchords);
	zs=repmat(mpx.bot.geom.slit.Z,1,def.bot.nchords);
	rd=mpx.bot.geom.wires.R;
	zd=mpx.bot.geom.wires.Z;
	phid=zeros(1,def.bot.nchords);
	tvd=zeros(1,def.bot.nchords);
	pvd=atan((zs-zd)./(rs-rd));
	pvd(pvd<0)=pi+pvd(pvd<0);
	pvd=pi-pvd;
	mpx.bot.geom.dcd=psitbxdcd(rd,zd,phid,pvd,tvd,nd);
    if (strcmp(detec,'bot'))
        rd=rd(chords-def.top.nchords);
        zd=zd(chords-def.top.nchords);
        phid=phid(chords-def.top.nchords);
        pvd=pvd(chords-def.top.nchords);
        tvd=tvd(chords-def.top.nchords);
    end
end
if(strcmp(detec,'both'))
	rs=[repmat(mpx.top.geom.slit.R,1,def.top.nchords) repmat(mpx.bot.geom.slit.R,1,def.bot.nchords)];
	zs=[repmat(mpx.top.geom.slit.Z,1,def.top.nchords) repmat(mpx.bot.geom.slit.Z,1,def.bot.nchords)];
    rd=[mpx.top.geom.wires.R mpx.bot.geom.wires.R];
	zd=[mpx.top.geom.wires.Z mpx.bot.geom.wires.Z];
	phid=zeros(1,def.top.nchords+def.bot.nchords);
	tvd=zeros(1,def.top.nchords+def.bot.nchords);
	pvd=atan((zs-zd)./(rs-rd));
	pvd(pvd<0)=pi+pvd(pvd<0);
	pvd=pi-pvd;
    % restrict to selected chords
    rd=rd(chords);
    zd=zd(chords);
    phid=phid(chords);
    pvd=pvd(chords);
    tvd=tvd(chords);
end
mpx.geom.dcd=psitbxdcd(rd,zd,phid,pvd,tvd,nd);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Load the DMPX filters       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'f'
    fprintf('\ncase f\n');
mdsopen(shot);
mdsdata('reset_public()');
fprintf('- Load the MPX filters and detection gas.\n');
if shot<26765,
	disp('No mobile absorber holder for shot<26765')
	mpx.top.filters.mobile_holder_pos=NaN;
	mpx.bot.filters.mobile_holder_pos=NaN;
	I_pos=4;
else
	if (exist('pos','var')==1 && sum(pos==[1 2 3 4])==1)
		disp(['Mobile absorber position forced to ' num2str(pos)])
		I_pos=pos;
	else
		disp('Load the mobile absorber holder position')
		a0=strcmp(deblank(mdsdata('\VSYSTEM::TCV_PUBLICDB_B["ILOT_NEB:A0"]')),'OFF');
		a1=strcmp(deblank(mdsdata('\VSYSTEM::TCV_PUBLICDB_B["ILOT_NEB:A1"]')),'OFF');
		a2=strcmp(deblank(mdsdata('\VSYSTEM::TCV_PUBLICDB_B["ILOT_NEB:A2"]')),'OFF');
		a3=strcmp(deblank(mdsdata('\VSYSTEM::TCV_PUBLICDB_B["ILOT_NEB:A3"]')),'OFF');
		I_pos=4*a0+3*a1+2*a2+1*a3;
	end
	mpx.top.filters.mobile_holder_pos=I_pos;
	mpx.bot.filters.mobile_holder_pos=I_pos;
	abso.material=abso.material{I_pos};
	abso.thickness=abso.thickness(I_pos);
end
I=find(shot>=[20030 21545 22043 23578 31942 32039]); %gas in the wire chamber
switch I(end) 
case 1 	       
	mpx.top.gas.type='ArCH';
	mpx.top.gas.thickness=6e-3;
case 2
	mpx.top.gas.type='KrCH';
	mpx.top.gas.thickness=6e-3;
case 3
	mpx.top.gas.type='XeCH';
	mpx.top.gas.thickness=6e-3;
case {4,6}
	if shot<26555
		mpx.top.gas.type='KrCH';
		mpx.top.gas.thickness=6e-3;
	else
		if strcmp(detec,'top')|strcmp(detec,'both'), 
			mpx.top.gas.type='KrCH';
			mpx.top.gas.thickness=8e-3;
		end
		if strcmp(detec,'bot')|strcmp(detec,'both'), 
			mpx.bot.gas.type='KrCH';
			mpx.bot.gas.thickness=7.6e-3; % bottom detector is slightly thiner (construction flaw...)
		end
	end  
case 5
	if strcmp(detec,'top')|strcmp(detec,'both'), 
		mpx.top.gas.type='ArCH';
		mpx.top.gas.thickness=8e-3;
	end
	if strcmp(detec,'bot')|strcmp(detec,'both'), 
		mpx.bot.gas.type='ArCH';
		mpx.bot.gas.thickness=7.6e-3; % bottom detector is slightly thiner (construction flaw...)
	end  
end
I=find(shot>=[20030 23806 26555]);
switch I(end) 
case 1
	mpx.top.filters.material={'He','Be'};
	mpx.top.filters.thickness=[526.8 350e-6];
	if shot>=21946 & shot<22085 %Al absorber inserted manually in the He tube
		mpx.top.filters.material{end+1}='Al';
		mpx.top.filters.thickness(end+1)=[308e-6];
	end 
case 2
	mpx.top.filters.material={'He','Be'};
	mpx.top.filters.thickness=[223.6e-3 150e-6];     
case 3
	if strcmp(detec,'top')|strcmp(detec,'both'), 
		mpx.top.filters.material={'He','Be'};
		mpx.top.filters.thickness=[229.6e-3 100e-6];
		%
		% Add the filter of the mobile absorber
		%
		if ~strcmp(abso.material,'nothing'),
			I=0;
			for iii=1:length(mpx.top.filters.material),
				if strcmp(abso.material,mpx.top.filters.material{iii})
					I=iii;
				end
			end
			if I>0,
				mpx.top.filters.thickness(I)=mpx.top.filters.thickness(I)+abso.thickness;
			else
				mpx.top.filters.material{end+1}=abso.material;
				mpx.top.filters.thickness(end+1)=abso.thickness;
			end
		end
	end
	if strcmp(detec,'bot')|strcmp(detec,'both'), 
		mpx.bot.filters.material={'He','Be','Air',mpx.top.gas.type};
		mpx.bot.filters.thickness=[229.6e-3 200e-6 10e-3 mpx.top.gas.thickness];
		%
		% Add the filter of the mobile absorber
		%
		if ~strcmp(abso.material,'nothing'),
			I=0;
			for iii=1:length(mpx.bot.filters.material),
				if strcmp(abso.material,mpx.bot.filters.material{iii})
					I=iii;
				end
			end
			if I>0,
				mpx.bot.filters.thickness(I)=mpx.bot.filters.thickness(I)+abso.thickness;
			else
				mpx.bot.filters.material{end+1}=abso.material;
				mpx.bot.filters.thickness(end+1)=abso.thickness;
			end
		end
		%
		% Add the additional absorber between the two detectors (if present)
		%
		I=find(shot>=[28751;28834;
				28862;28891;
				28943;28973
				30681;30690]); % Insert shots by pair, additional absorber betweem the two shots numbers
		addi_material={'Al','Al','Al','Al'};
		addi_thick=[308 308 308]*1e-6;
		if ~isempty(I)&rem(I(end),2)~=0, % Means there was an additional absorber between the 2 detectors
			II=0;
			for iii=1:length(mpx.bot.filters.material),
				if strcmp(addi_material(rem(I(end),2)),mpx.bot.filters.material{iii}),
					II=iii;
				end
			end
			if II>0,
				mpx.bot.filters.thickness(II)=mpx.bot.filters.thickness(II)+addi_thick(rem(I(end),2));
			else
				mpx.bot.filters.material{end+1}=addi_material{rem(I(end),2)};
				mpx.bot.filters.thickness(end+1)=addi_thick(rem(I(end),2));
			end
		end
	end	 
end
mdsclose;
% fprintf('- MDS disconnect.\n');
% mdsdisconnect;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Load the detectors efficiency       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'e'
    fprintf('\ncase e\n');
if strcmp(detec,'top')|strcmp(detec,'both'),
	fprintf('- Load the top detector efficiency.\n');
	s='';
	for iii=1:length(mpx.top.filters.material)
		s=[s mpx.top.filters.material{iii} num2str(mpx.top.filters.thickness(iii)*1e6) ', '];
	end
	s(end-1:end)=' $';
	s=[s ' ' mpx.top.gas.type 'abs' num2str(mpx.top.gas.thickness*1e6)];	 
	c_dir=pwd;
	filelocate=which('mpxdata.m');
	ind=findstr(filelocate,'mpxdata.m');
	if or(strcmp(computer,'IBM_RS'),strcmp(computer,'GLNXA64'))
		directory=[filelocate(1:ind-1) 'mpx_efficiency/'];
	else
		directory=[filelocate(1:ind-1) 'mpx_efficiency\'];
	end     	  
	eval(['cd ' directory]);
	out=XraySignals(s);
	out.response(isnan(out.response))=0;
	eval(['cd ' c_dir]);
	mpx.top.eff.ev=out.ev;
	mpx.top.eff.response=out.response;
end
if strcmp(detec,'bot')|strcmp(detec,'both'),
	fprintf('- Load the bottom detector efficiency.\n');
	s='';
	for iii=1:length(mpx.bot.filters.material)
		s=[s mpx.bot.filters.material{iii} num2str(mpx.bot.filters.thickness(iii)*1e6) ', '];
	end
	s(end-1:end)=' $';
	s=[s ' ' mpx.bot.gas.type 'abs' num2str(mpx.bot.gas.thickness*1e6)]; 
	c_dir=pwd;
	filelocate=which('mpxdata.m');
	ind=findstr(filelocate,'mpxdata.m');
	if or(strcmp(computer,'IBM_RS'),strcmp(computer,'GLNXA64'))
		directory=[filelocate(1:ind-1) 'mpx_efficiency/'];
	else
		directory=[filelocate(1:ind-1) 'mpx_efficiency\'];
	end     	  
	eval(['cd ' directory]);
	out=XraySignals(s);
	out.response(isnan(out.response))=0;
	eval(['cd ' c_dir]);
	mpx.bot.eff.ev=out.ev;
	mpx.bot.eff.response=out.response;	 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Load rho for each DMPX chord       %	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'r'
    fprintf('\ncase r\n');
fprintf('- Calculate rho_psi for each detector chord.\n');
if(~exist('opt','var'))
	opt=def.opt;
end
if(exist('psi','var')~=1)
    if ~exist('rtime','var'), 
        if(exist('time','var') && ~ischar(time))
            rtime=time;
        else
            rtime=def.rtime;
        end
    end
    if(~exist('liuqe','var')) % problematic: || isempty(find(repmat(liuqe,1,3)==[1 2 3])))
		liuqe=def.liuqe;
    elseif(~any(liuqe==[1 2 3]))
		liuqe=def.liuqe;
    end
	mdsopen(shot);
    mdsdata('reset_public()');
    psi_t=mdsdata('dim_of(\results::psitbx:psimag)'); %to get the psitbx time base
	mdsclose;
%     fprintf('- MDS disconnect.\n');
%     mdsdisconnect;
    if(isempty(psi_t))
        warning('No psitbx data');
        continue; % Go to next action.
    end
	I=iround(psi_t,rtime);
	t=psi_t(I(1):I(end));
	disp('Loading poloidal flux with PSITBXTCV...');
    if(liuqe==1)
		psi = psitbxtcv(shot,t, '01');
 	else
		psi = psitbxtcv(shot,t, '01', ['LIUQE' num2str(liuqe)]);   
    end
else
	t=psi.t; 	 
end 
if (exist('vol','var') && vol~=0) 
	fsd = psitbxp2p(psi,'FS'); 
	fsg = psitbxfsg(fsd);
	s=size(fsg.vol.x);	  
	rhovol=sqrt(fsg.vol.x./repmat(max(fsg.vol.x),s(1),1));
else
	vol=def.vol;
end
if strcmp(detec,'top')|strcmp(detec,'both'), 
	mpx.top.rho.rhopsi=rhochords(mpx.top.geom.dcd,psi);
	mpx.top.rho.time=t;
	if vol~=0,
		mpx.top.rho.rhovol=NaN.*zeros(length(t),def.top.nchords);
		for iii=1:length(t),
			%mpx.top.rho.rhovol(ii,:)=interpos(3,fsg.vol.grid.x{1},rhovol(:,ii),mpx.top.rho.rhopsi(ii,:));
            mpx.top.rho.rhovol(iii,:)=spline(fsg.vol.grid.x{1},rhovol(:,iii),mpx.top.rho.rhopsi(iii,:));
		end
	end
	if opt<0 % rho from -1 to 1
		sr=repmat(mpx.top.geom.slit.R,1,length(t));	   
		sz=repmat(mpx.top.geom.slit.Z,1,length(t));
		theta_mag=atan((psi.zmag-sz)./(psi.rmag-sr));
		%if theta_mag<0, theta_mag=theta_mag+pi; end; % Doesn't work !
		theta_mag(theta_mag<0)=theta_mag(theta_mag<0)+pi; % Correction L. Curchod (15.07.2008)
		theta_mag=theta_mag-pi/2; %angle between the vertical and the line between the MPX slit center and the magnetic axis, trigonometric direction (=negative values when R_axis>0.88)
		I_HFS=repmat(pi/2-mpx.top.geom.dcd.pvd,length(t),1)>repmat(theta_mag',1,length(mpx.top.geom.wires.R)); %select the HFS chords
		mpx.top.rho.rhopsi(I_HFS)=-mpx.top.rho.rhopsi(I_HFS);
		if vol~=0
			mpx.top.rho.rhovol(I_HFS)=-mpx.top.rho.rhovol(I_HFS);
		end
	end
   
end	 
if strcmp(detec,'bot')|strcmp(detec,'both'), 
	mpx.bot.rho.rhopsi=rhochords(mpx.bot.geom.dcd,psi);
	mpx.bot.rho.time=t;
	if vol~=0,
		mpx.bot.rho.rhovol=NaN.*zeros(length(t),def.bot.nchords);
		for iii=1:length(t),
			%mpx.bot.rho.rhovol(ii,:)=interpos(3,fsg.vol.grid.x{1},rhovol(:,ii),mpx.bot.rho.rhopsi(ii,:));
            mpx.bot.rho.rhovol(iii,:)=spline(fsg.vol.grid.x{1},rhovol(:,iii),mpx.bot.rho.rhopsi(iii,:));
		end
	end
	if opt<0 % rho from -1 to 1 
		sr=repmat(mpx.bot.geom.slit.R,1,length(t));	   
		sz=repmat(mpx.bot.geom.slit.Z,1,length(t));
		theta_mag=atan((psi.zmag-sz)./(psi.rmag-sr));
		%if theta_mag<0, theta_mag=theta_mag+pi; end; % Doesn't work !
		theta_mag(theta_mag<0)=theta_mag(theta_mag<0)+pi; % Correction L. Curchod (15.07.2008)
		theta_mag=theta_mag-pi/2; %angle between the vertical and the line between the MPX slit center and the magnetic axis, trigonometric direction (=negative values when R_axis>0.88)
		I_HFS=repmat(pi/2-mpx.bot.geom.dcd.pvd,length(t),1)>repmat(theta_mag',1,length(mpx.bot.geom.wires.R)); %select the HFS chords
		mpx.bot.rho.rhopsi(I_HFS)=-mpx.bot.rho.rhopsi(I_HFS);
		if vol~=0
			mpx.bot.rho.rhovol(I_HFS)=-mpx.bot.rho.rhovol(I_HFS);
		end
	end	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Unknown action       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

otherwise
error('mpxdata:UnknownAction',['''' action(ii) '''' ' is an unknown action'])
end % End of 'for' loop
end % End of 'switch' cases
%disp(['Closing MDS connection to shot ' num2str(shot) '.']);
mdsclose;
%mdsdisconnect;

%path(backup_path);

fprintf('\n+ MPXDATA finished\n\n');

end