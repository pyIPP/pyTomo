% find_dR is used for the minimisation of the difference between
%	the profile obtained by the HFS chords and the one obtained by
%	the LFS chords.
%
%
%[param]=find_dR(dR,fsdmod,psi,I_chords_HFS,I_chords_LFS,MPX_signal,MPX_sig_err,geom,meth)
%
%
%See also
%
%	interpos.m  prof_inv.m
%
%
	
function [param]=find_dR(dR,fsdmod,psi,I_chords_HFS,I_chords_LFS,MPX_signal,MPX_sig_err,geom,meth)

%normalisation of parameter dR to help fminsearch convergence
dR=dR/500;
dR(1)=dR(1)/5;

%edge conditions check
if max(abs(dR)) > 0.04,
	param = realmax;
	return
end

%disp([num2str(dR(1)) '   ' num2str(dR(2))])
%inversion for each type of chords
[rho_HFS,g_HFS,chi2_HFS,T_HFS]=prof_inv(psi,[fsdmod dR],I_chords_HFS,MPX_signal,MPX_sig_err,geom,meth);
[rho_LFS,g_LFS,chi2_LFS,T_LFS]=prof_inv(psi,[fsdmod dR],I_chords_LFS,MPX_signal,MPX_sig_err,geom,meth);

if any(isnan(g_HFS))|any(isnan(g_LFS)),
	param = realmax;
	return
end

%interpolate on the same rho basis
rhotcv=[0:0.025:1];
gg_HFS=interpos(13,rho_HFS,g_HFS,rhotcv); %line vector
gg_LFS=interpos(13,rho_LFS,g_LFS,rhotcv);


tmp1=gg_HFS-gg_LFS;
tmp2=mean([gg_HFS;gg_LFS]);
tmp2(tmp2<max(tmp2)/10)=max(tmp2)/10; %diminue le poids du bord 
tmp=(tmp1./tmp2).^2;
param = log(0.01+sum(tmp)); %help fminsearch by increasing the gradients near the minimum 
