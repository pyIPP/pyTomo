%plas_grid.m returns a matrix of ones and zeros based on the dimension of
%	the grid contained in psi. The value is one for a point of the grid 
%	contained in the last closed flux surface.
%
%
%A = plas_grid(psi)
%
%
%Inputs:
%
%	psi	psitbxpsi object of format '01' at a unique time.
%
%
%Outputs:
%
%	A	matrix of 1 (inside the plasma), 0 (outside)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Jonathan Rossel, 17.02.05


function A = plas_grid(psi)
	
if psi.format ~= '01',
	disp('psi does not have the good format')
	return
end
if length(psi.psitbxfun.t) ~= 1,
	disp('the equilibrium must be found for one time only')
end

%flux surface description object
fsd = psitbxp2p(psi,'FS');

%grid of psi in flux coordinates
PsiOnFlux = psitbxg2g(psi.psitbxfun.grid,'Flux',psi);

%selects the points situated inside the last closed flux surface
tmp = find(PsiOnFlux.x{1}<=1);
PsiOnFlux = psitbxgrid('Flux','Points',{PsiOnFlux.x{1}(tmp),PsiOnFlux.x{2}(tmp)},psi);

%return in cylindrical coordinates with a subset of points and finds the points of the original grid
SmallCylGrid = psitbxg2g(PsiOnFlux,'Cylinder',fsd);
tmp = find(~isnan(SmallCylGrid.x{1}));
tmpR = SmallCylGrid.x{1}(tmp);
tmpZ = SmallCylGrid.x{2}(tmp);
tmp = iround(psi.psitbxfun.grid.x{1},tmpR,1);
tmp2 = iround(psi.psitbxfun.grid.x{2},tmpZ,1);
tmpR = psi.psitbxfun.grid.x{1}(tmp)';
tmpZ = psi.psitbxfun.grid.x{2}(tmp2)';
tmp = sortrows([tmpR,tmpZ]);
In_points = unique(tmp,'rows');
[tmpR,tmpZ] = meshgrid(psi.psitbxfun.grid.x{1},psi.psitbxfun.grid.x{2});
Points = [reshape(tmpR,size(tmpR,1)*size(tmpR,2),1),reshape(tmpZ,size(tmpZ,1)*size(tmpZ,2),1)];
[tmp,ind,tmp] = intersect(Points,In_points,'rows');
tmp = zeros(size(tmpR));
tmp(ind) = 1;
tmp = tmp';
A = zeros(size(psi));
A(psi.psitbxfun.x<=1) = 1;
A = A.*tmp;
