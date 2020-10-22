%Calculate the rho_psi for each diagnostic chord
%
% function [rho_psi]=rho_chords(chordR,chordZ,psi)
%
% Output :	rho_psi -> rho minimum for each chord
%
% Input :	chordR -> R-coordinates in m (two points for each chord), dim : 2*nb_chords
% 		chordZ -> Z-coordinates in m (two points for each chord), dim : 2*nb_chords
% 		psi -> psi=psitbxtcv(shot,time,'01') poloidal flux in cylindrical coordinates
%		dRZ_mag -> (optional) modification of the position of the magnetic axis :
%			dRZ_mag(1)=R_new-R_liuqe
%			dRZ_mag(2)=Z_new-Z_liuqe


function [rho_psi]=rho_chord(chordR,chordZ,psi,dRZ_mag)
global_p;

if nargin<4
 dRZ_mag=[0 0];
end;


%MPX chords description in a psitbxdcd object
rd=chordR(1,:)-dRZ_mag(1);	%il serait plus elegant de changer r_mag et z_mag dans 'psi', mais l'acces a la structure est protege...
zd=chordZ(1,:)-dRZ_mag(2);
phid=zeros(1,64);
tvd=zeros(1,64);
nd=1000;
pvd=atan((chordZ(2,:)-chordZ(1,:))./(chordR(2,:)-chordR(1,:)));
pvd(pvd<0)=pi+pvd(pvd<0);
pvd=pi-pvd;

dcd=psitbxdcd(rd,zd,phid,pvd,tvd,nd);
GridMPXChords=psitbxgrid('Diagnostic-Chords','Points',{intersect(dcd,psi.grid)},dcd);
GridMPXChordsFlux = psitbxg2g(GridMPXChords,'Flux',psi)

rho_psi=min(reshape(GridMPXChordsFlux.x{1},nd,64),[],1);
rho_psi(rho_psi>1)=NaN;


%Calcule le rho_psi (point de tangence) pour chaque corde. Renvoie NaN si la corde est en dehors du plasma.
%%Methode 1
%n=1000
%N=64
%R=zeros(n,N);
%Z=zeros(n,N);
%for ii=1:N,
%  R(:,ii)=linspace(chordR(1,ii),chordR(2,ii),n)';
%  Z(:,ii)=linspace(chordZ(1,ii),chordZ(2,ii),n)';
%end;
%chords=psitbxgrid('Cylindrical','Points',{[R] [Z] [0]}) %description of diagnostics chords in cylindrical coordinates
%rho_chords=psitbxg2g(chords,'Flux',psi); %description of diagnostics chords in Flux-coordinates
%rho_chords=min(rho_chords.x{1},[],1);

