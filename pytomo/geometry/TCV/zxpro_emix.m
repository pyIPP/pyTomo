function emix=zxpro_emix(Z,te,thick)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function emix=emix(Z,te,Bethick)
%   supplies normalized emissivity for z and given te vector
%   preferred use with global definition in main program to avoid 
%   frequent disk access:
%	  global ION_Z POLYORDER_0 POLYORDER_47 IONPOLY_0 IONPOLY_47
% 	 load ionmatrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==2,thick=47;end % default Be foil for tomography
scale=4e-32; % scaling from IONEQ defaults
global ION_Z POLYORDER_0 POLYORDER_47 IONPOLY_0 IONPOLY_47 

if length(ION_Z)==0, load zxpro_ionmatrix,end % casual use, normally not necessary
index=find(ION_Z==Z);
if thick~=47 & thick~=0,
 disp(['EMIX error: requested filter thickness unavailable', int2str(thick),' available:0, 47']) 
end
ok=find(te>=200 & te<=20000 &te~=NaN);
emix=NaN*ones(size(te));
if thick==47
 p=IONPOLY_47(index,15-POLYORDER_47(index):15);
else
 p=IONPOLY_0(index,15-POLYORDER_0(index):15);
end
   emix(ok)=scale*exp(polyval(p,log(te(ok))));
 
