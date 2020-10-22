% Calculation of the local emissivity profile from line integrated measurement,
% assuming constant emissivity on a flux surface. There are several ways to
% use this routine:
%
% (1) The flux surface description (FSD) modification is specified: the
% inversion is done for each time of the signal (keeping a constant transfer matrix).
%
% (2) The flux surface description (FSD) modification is NOT specified and
% size(signal.dim{1}) = 1, that is only 1 time (typically, signal.data is the
% signal manually averaged on the time interval of interest): the routine computes the
% parameters for the FSD modification so as to minimize the difference between
% the profiles obtained with the HFS and the LFS chords. Then it stops. You can check
% the validity of the FSD modification and then use the resulting FSD modification
% parameters for the inversion on each time of the interval of interest (go to
% case (1)).
%
% (3) The flux surface desciption (FSD) modification is NOT specified and
% size(signal.dim{1}) > 1, that is more than 1 time: the signal is first averaged
% on the whole time interval and the FSD modification parameters are computed so as
% to minimize the difference between the profiles obtained with the HFS and the
% LFS chords of the DMPX. Then the routine DOESN'T stop and the inversion is done
% for each time of the non-averaged signal. It is thus faster but you cannot check
% the result of the FSD modification.
%
% It is advised to use case (2): first compute the FSD modification parameters for
% the signal AVERAGED on a time interval during which the plasma is stationnary (position,
% density, B field, plasma current,...) or over several MHD activity periods. Then, use
% these parameters to do the inversion for each time.
%
% [fsd_mod,new_psi,rho,g,chi2,T,I_chords,rho_HFS,g_HFS,chi2_HFS,T_HFS,I_chords_HFS,rho_LFS, ...
% 
% g_LFS,chi2_LFS,T_LFS,I_chords_LFS]=mpx_inv(signal,sig_err,geom,psi,ok_chords,fsd_mod,inv_method)
%
% Inputs:
%	signal: must have two fields.
%		signal.dim{1}: time
%		signal.dim{2}: chords number
%		signal.data: (time*chords)
%	sig_err: relative error on signal (same dimension as ok_chords)
%	geom: detector geometry, as obtained by mpxdata(shot,'g'); 
%	psi: flux surface description, as obtained by psi = psitbxtcv(shot,time,'01');
%	     (must have length(time)=1)
%	ok_chords: array with 1 if the chord is ok, 0 if not
%		   (default: ok_chords=ones(1,length(signal.dim{2}));)
%	fsd_mod: parameters for the flux surface description modification (default: fsd_mod=3;)
%		0	no modification 
%		[1,dR]	radial shift of the whole equilibrium (dR>0, shift in the LFS direction)
%		[2,dR]	radial shift proportional to dR*(1-rhopsi) (deformation of the equilibrium)
%		[3,dR1,dR2] radial shift of dR1 + dR2*(1-rhopsi)
%		[4,dR,dZ] radial and vertical shift
%		If length(fsd_mod)=1, the minimization process is used to find the flux surface
%		description modification.
%		If fsd_mod(1)<0, the minimization process is launched using fsd_mod(2:end) as
%		first guess parameters for dR.
%	inv_method: inversion method used
%		1:	mixes the minimum Fisher information and the
%			minimum second order derivative 
%		2: 	minimum second order derivative
%		3:	minimum Fisher information (default - best one)
%		4:	least square with Npixel = Nchords/2
%
% Outputs:
%	fsd_mod: parameters for the flux surface description modification (default: fsd_mod=3;);
%		 returns the input fsd_mod, if given.
%	new_psi: new flux surface description, modified with fsd_mod
%	rho: rho_psi of inverted profiles
%	g: inverted profiles
%	chi2: chi2=mean((T*g-f).^2) that is the mean square difference between the resulting
%	      inverted profiles and the initial integrated data, at each time.
%	      !!! IF CHI2 > 1, THE RESULT IS COMPLETELY FALSE !!!
%	T: transfer matrix
%	I_chords: chords actually used for the inversion (not necesseraly equal to ok_chords)
%	
%	The following outputs whith a name ending with _HFS (resp. _LFS) are the same outputs
%	than above but for the inverted profile of averaged data from the HFS (resp. LFS) chords only.
%	


function [fsd_mod,new_psi,rho,g,chi2,T,I_chords,rho_HFS,g_HFS,chi2_HFS,T_HFS,I_chords_HFS,rho_LFS,g_LFS,chi2_LFS,T_LFS,I_chords_LFS]=...
	mpx_inv(signal,sig_err,geom,psi,ok_chords,fsd_mod,inv_method)

addpath /home/matlab/crpptbx-7.6.0/dmpx/mpx_inversion

%%%%%%%%%%%%% Inputs check %%%%%%%%%%%%%
if nargin<4,
 error('mpx_inv:WrongInput',['Inputs ''signal'', ''sig_err'', ''geom'' and ''psi'' required'])
end
if ~(isfield(signal,'dim')&isfield(signal,'data'))
 error('mpx_inv:WrongInput',['Input ''signal'' has not the required fields'])
end
if length(psi.t)~=1, 
 error('mpx_inv:WrongInput',['Input ''psi'' must be given for only one time'])
end
if exist('ok_chords')~=1, %default value for ok_chords.
 ok_chords=ones(1,length(signal.dim{2}));
elseif isempty(find(length(ok_chords)==size(signal.data)))
 error('remove_spike:WrongInput',['Input ''signal'' and ''ok_chords'' must have one dimension in common']);
end
if exist('fsd_mod')~=1, %default value for fsd_mod.
 fsd_mod=3;
end
minimize=0;
if fsd_mod(1)<0,
	fsd_mod(1)=abs(fsd_mod(1));
	minimize=1; %then minimization will be done with fsd_mod(2:end) as first guess parameters for dR.
end
if length(fsd_mod)~=1,
 switch fsd_mod(1),
  case 1|2
   if length(fsd_mod)~=2, error('mpx_inv:WrongInput',['Input ''fsd_mod'' has not the required length']); end;
  case 3|4
   if length(fsd_mod)~=3, error('mpx_inv:WrongInput',['Input ''fsd_mod'' has not the required length']); end;
 end
else
 if fsd_mod(1)~=0, minimize=1; end %then minimization will be done.
end
if  exist('inv_method')~=1, %default value for inv_method.
 inv_method=3;
end


%%%%% MPX signal normalisation %%%%%
A_fente=geom.slit.dR.*geom.slit.length;
A_fil=(geom.wires.R(1)-geom.wires.R(2)).*geom.wires.length; %Aire efficace de chaque fil = aire de la zone de detection au niveau du fil. Efficient area of each wire = area of the detection zone on the wire level.

geo_fact=A_fente.*A_fil.*cos(pi/2-geom.dcd.pvd).^2./(4.*pi.*(geom.slit.Z-geom.wires.Z).^2);    
%Rajouter l'aire de la fente pour chaque fil (obstruction des tuiles pour certains fils) -->Pochon, Chavan, inclure dans mpxdata.m

signal.data=signal.data./repmat(geo_fact,size(signal.dim{1},1),1);



%%%%%%%%%%%%% Let's go %%%%%%%%%%%%%
if ~any(size(signal.data)==1), 
	MPX_signal=mean(signal.data);
else
	MPX_signal=signal.data;
end
I_chords=[1:length(signal.dim{2})]; 	
I_chords_ref=[1:length(signal.dim{2})];
I_chords(isnan(MPX_signal)|~ok_chords)=[];
disp(['Chords removed: ' num2str(find(isnan(MPX_signal)|~ok_chords))])


%initialisation for dR
if minimize==1,	
	if length(fsd_mod)==1,		%first guess of dR using maximum signal in DMPX profile
		[tt,ii]=max(interpos(13,I_chords,MPX_signal(I_chords),I_chords,1));
		dR=geom.slit.R + (psi.zmag-geom.slit.Z)*tan(geom.dcd.pvd(I_chords(ii))-pi/2) - psi.rmag;
		if abs(dR) > 0.05, dR = 0; end
		switch fsd_mod(1) 
			case 1
				fsd_mod(2)=dR;
			case 2
				fsd_mod(2)=dR;
			case 3
			%mdsopen(shot)
 			%tmp_time=mdsdata('\results::xtoma:time');
			%if ~isempty(tmp_time)
			%	tmp=mdsdata('\results::xtoma:rcat');
			%	I=iround(tmp_time,[signal.dim{1}(1) signal.dim{1}(end)])
 			%	tmp=mean(tmp(I(1):I(2)))/100;
			%	dR=tmp-psi.rmag;
			%	fsd_mod(2:3)=[0 dR];
			%else
				fsd_mod(2:3)=[dR/2 dR/2];
			%end
			case 4
				fsd_mod(2:3)= [dR 0];
			otherwise
				error('mpx_inv:WrongInput',['Input ''fsd_mod(1)'' has an unauthorized value'])
		end
	end
end


%flux surface modification
if fsd_mod(1)~=0,
	res=fsd_mod(2:end);
	res(res==0)=5e-4; %if dR guess is set to 0, fminsearch doesn't make big steps in dR. Strange...
	res = res*500;  %rescale input parameters for fminsearch, increases the convergence
	res(1)=res(1)*5; %moves 5 times faster in dR2 direction than dR1 in fminsearch
	options = optimset('MaxIter',80,'TolX',0.01*500,'TolFun',0.001,'Display','iter');
	new_iter=1;
	i_count=0;
	disp('please wait...')
	while new_iter==1
		i_count=i_count+1;
		new_psi = psi_modif(psi,fsd_mod);  
		theta_mag = atan((new_psi.zmag-geom.slit.Z)/(new_psi.rmag-geom.slit.R));
		if theta_mag<0
			theta_mag=theta_mag+pi;
		end
		theta_mag=pi-theta_mag; %standard psitbx angle def (0 at HFS and clock-wise orientation)
	
		I_chords_HFS=I_chords_ref((geom.dcd.pvd<theta_mag)&~isnan(MPX_signal)&ok_chords);
		I_chords_LFS=I_chords_ref((geom.dcd.pvd>=theta_mag)&~isnan(MPX_signal)&ok_chords);
	
		if minimize==1,		
			res=fminsearch(@find_dR,[res],options,fsd_mod(1)+10,psi,I_chords_HFS,I_chords_LFS,MPX_signal,sig_err,geom,inv_method);
			options=optimset('MaxIter',80,'TolX',0.002*500,'TolFun',0.0005,'Display','iter');
			fsd_mod(2:end)=res/500;
			fsd_mod(2)=fsd_mod(2)/5;
			if i_count==3
				new_iter=0;
			end
		else
			new_iter=0;
		end
	
	end
end

new_psi = psi_modif(psi,fsd_mod);  
[rho,g,chi2,T,I_chords]=prof_inv(psi,fsd_mod,I_chords,signal.data,sig_err,geom,inv_method);
[rho_HFS,g_HFS,chi2_HFS,T_HFS,I_chords_HFS]=deal(NaN);
[rho_LFS,g_LFS,chi2_LFS,T_LFS,I_chords_LFS]=deal(NaN);

% pour pouvoir controler la qualite de la modification des surfaces de flux
	theta_mag = atan((new_psi.zmag-geom.slit.Z)/(new_psi.rmag-geom.slit.R));
	if theta_mag<0, theta_mag=theta_mag+pi; end;
	theta_mag=pi-theta_mag; %standard psitbx angle def (0 at HFS and clock-wise orientation)
	I_chords_HFS=I_chords_ref((geom.dcd.pvd<theta_mag)&~isnan(MPX_signal)&ok_chords);
	I_chords_LFS=I_chords_ref((geom.dcd.pvd>=theta_mag)&~isnan(MPX_signal)&ok_chords);

	if length(I_chords_HFS)>2,
		[rho_HFS,g_HFS,chi2_HFS,T_HFS,I_chords_HFS]=prof_inv(psi,fsd_mod,I_chords_HFS,MPX_signal,sig_err,geom,inv_method);
	end
	if length(I_chords_LFS)>2,
		[rho_LFS,g_LFS,chi2_LFS,T_LFS,I_chords_LFS]=prof_inv(psi,fsd_mod,I_chords_LFS,MPX_signal,sig_err,geom,inv_method);
	end
