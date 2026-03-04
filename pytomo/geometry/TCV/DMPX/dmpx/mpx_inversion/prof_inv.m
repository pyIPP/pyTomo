%
%[rho,g,chi2,T,I_chords]=prof_inv(psi,fsd_mod,I_chords,signal,sig_err,geom,meth);
%
%
%Inputs:
%
%	psi:		psitbxpsi object, format '01', a unique time
%	fsd_mod:		parameter of modification of the poloidal flux (see mpx_inv.m)
%	I_chords:		chords of interest
%	signal:		signal (dim: time*chords)
%	sig_err:		absolute error on the signals
%	geom:	structure containing the geometry of the DMPX
%	meth:		inversion method (see mpx_inv_stat.m)
%
%
%Outputs:
%
%	rho		radial flux coordinate of the pixel centers in rho_psi
%	g		solution emissivity g(rho)
%	chi2		error parameter given by the inversion methods
%	T		transfer matrix (see mpx_transf_mat.m)
%	I_chords	chords of the input I_chords that pass inside the plasma
%
%


function [rho,g,chi2,T,I_chords]=prof_inv(psi,fsd_mod,I_chords,signal,sig_err,geom,meth)


N_chords=length(I_chords);


%check that the selected chords are in the plasma
rho_FS = [0 1];
psi_new = psi_modif(psi,fsd_mod); %changes the flux description according to fsd_mod
[sl,sh,pc] = psitbxdxp(geom.dcd,psi_new,rho_FS); %calculate the intersection of the chords with the flux surface
						%sl and sh contains the chords-coordinates of the intersection (low and high)
						%pc contains the value of the psi01 flux at the closest point of the chord to the magnetic axis. Nan if not in the plasma
						%psi_modif is always needed to remove the problematic of the flux in the divertor zone.


I_ok = zeros(N_chords,1);
I_ok(~isnan(pc.x(I_chords))) = 1;

if  sum(~I_ok)>0
	%disp(['Chord(s) ' num2str(I_chords(~I_ok)) ' not in the plasma.']);
	I_chords=I_chords(I_ok==1);
	N_chords=length(I_chords);
end


%compute the transfer matrix
tmp = sort(sqrt(max(pc.x(I_chords),0))); % closest position of the chord in flux coord
tmp = tmp + [diff(tmp)/2,0];
rho_FS = [0,tmp(1:end-1),1];
rho_FS = unique(rho_FS); %exactly one pxl per chords, or less (could be used for direct inversion)
if meth ~= 4,
	rho_FS = [0,linspace(rho_FS(2),1,N_chords+1)]; %to assure that the central pixel has a chord inside, use +1 to be sure that there will be a regularisation with minfisher
else
	rho_FS = [rho_FS(1:2:end-1),1]; %takes one on two plus the end;
end
[sl,sh,pc] = psitbxdxp(geom.dcd,psi_new,rho_FS);
sl = squeeze(sl.x{1});
sh = squeeze(sh.x{1});

[T]=mpx_transf_mat(sl,sh,I_chords);


%inversion of the transfer matrix
switch meth
 case 1
  meth_name='FishSecMix';
 case 2
  meth_name='SndOrdReg';
 case 3
  meth_name='new_minfisher';
 case 4
  meth_name='min_chi2';
end

rho=(rho_FS(2:end)+rho_FS(1:end-1))/2;
n_t=size(signal,1);
g=zeros(length(rho),n_t);
chi2=zeros(1,n_t);
rho_out=zeros(length(rho),n_t);
sig_err_abs=signal.*repmat(sig_err,n_t,1);

l=0;
for ii=1:n_t,
	eval(['[g(:,ii),chi2(1,ii),tmp]=' meth_name '(signal(ii,I_chords)'',sig_err_abs(ii,I_chords)'',T,rho);']);
%	eval(['[g(:,ii),chi2(1,ii),tmp]=' meth_name '(signal(I_chords,ii),sig_err_abs(I_chords,ii),T,rho,0.01);']);
	rho_out(:,ii)=tmp';
	if rem(ii,10)==0, 
		for jj=1:l, fprintf(1,'\b'); end
		l=fprintf(1,'%3.1f%% done',ii/n_t*100);
	end;
end
if n_t>1, fprintf(1,'\n'); end
