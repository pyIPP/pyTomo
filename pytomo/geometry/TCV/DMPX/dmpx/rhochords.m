%Get rho_psi for each diagnostic chord using the psitoolbox (for a non moving detector)
%
% function [rho]=rhochords(dcd,psi)
%
%
% Input :	dcd:	diagnostic chord descritpion (psitbxdcd object)
%			psi:	flux description in cylindrical coordinates (psitbxpsi object)
%
% Output :	rho:	min(rho) obtained for each chord
%			(positive values, NaN if the chord is outside of the plasma)
%
% Remark: Here rho_psi=sqrt(1-psi/psi_max), as in the psitbx, GUIprofs,... but not as in LIUQE where rho_psi=1-psi/psi_max
%
% Example for the MPX: (will give rho from HFS to LFS, following the ordering of the channels in the dcd object)
%
%	mpx=mpxdata(shot,'g','detec','top');
%	dcd=mpx.top.geom.dcd;
%	psi=psitbxtcv(shot,time,'01');
%	rho=rhochords(dcd,psi)	


function [rho]=rhochords(dcd,psi)

GridMPXChords=psitbxgrid('Diagnostic-Chords','Points',{intersect(dcd,psi.grid)},dcd);
GridMPXChordsFlux = psitbxg2g(GridMPXChords,'Flux',psi);

rho=reshape(min(GridMPXChordsFlux.x{1},[],1),length(dcd.zd),length(psi.t))';
rho(rho>1)=NaN;


