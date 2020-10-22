function rho = get_mpx_rho(shot,time);
global_p;

load mpx_geometry.mat

psi = psitbxtcv(shot,time);

psi = psitbxp2p(psi,'01');      


[rd,zd,phid,pvd,tvd] = psitbx_xtomo_geometry(chordR/100,chordZ/100);


dcd = psitbxdcd(rd,zd,phid,pvd,tvd);

GridDiagChords = psitbxgrid('Diagnostic-Chords','Points', ...
		 {intersect(dcd,psi.grid)},dcd);
		 
GridDiagChordsFlux = psitbxg2g(GridDiagChords,'Flux',psi);


rho = min(GridDiagChordsFlux.x{1});

for ii=1:64
   r(1,ii,1)=rho(1,1,ii);
end

clear rho

%kkk1   = find(r>1);
%r(kkk1)= [];
kkk    = find(diff(r)<0);
lll    = find(diff(kkk)>1);

if isempty(lll)
   r(kkk)=-r(kkk);
else
   kkk(lll+1:end) = [];
   r(kkk)=-r(kkk);
end

%r(kkk)=-r(kkk);


rho = r;

