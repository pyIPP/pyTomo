%get the one time slice PsiTbx-Poloidal-Flux object (new_psi) corresponding to one time (t) from a multi-time slice PsiTbx-Poloidal-Flux object (psi)

%new_psi=onetime_psi(psi,t)

function new_psi=onetime_psi(psi,t)

if nargin<2,
	disp('Not enough input arguments')
	return
end
if length(t)>1,
	disp('Wrong input: length(t)~=1')
	return
end

new_psi=psi;
if length(psi.psitbxfun.t)>1,
	I=iround(psi.psitbxfun.t,t);
	new_psi=psitbxpsi(squeeze(psi.psitbxfun.x(:,:,I)),psi.psitbxfun.grid,psi.psitbxfun.t(I),'01');
	new_psi.rmag=psi.rmag(I);
	new_psi.zmag=psi.zmag(I);
	new_psi.psimag=psi.psimag(I);
end
