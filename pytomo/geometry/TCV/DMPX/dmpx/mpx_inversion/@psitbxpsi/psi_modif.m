%psi_modif.m modifies the poloidal flux (format '01'). The first modification
%	concerns the divertor zone: psi is set to 2-psi there. The other
%	modifications are relative to mpx_inv_stat.m.
%
%
%new_psi = psi_modif(psi,mod_p)
%
%
%Inputs:
%
%	psi		psitbxpsi object of format '01', at a given time.
%
%	mod_p		parameter of modification of psi:
%
%				0/10		only the divertor zone changes
%
%				[1/11, dR]	radial translation of dR
%
%				[2, dR]		radial deformation of dR. The
%						flux surfaces are translated by:
%						(1-rho)*dR, rho: flux coord. of
%						the surface. dR is tested to 
%						stay inside the plasma.
%
%				[12, dR]	idem but no test on dR
%
%				[3, dR1, dR2]	radial transl. of dR1, and
%						radial deform. of dR2.
%
%				[13, dR1, dR2]	idem but no test on dR2
%
%				[4/14, dR, dZ]	radial and vertical translation.
%
%
%Outputs:
%
%	new_psi		modified psitbxpsi object
%
%
%See also:
%
%	plas_grid.m
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Jonathan Rossel, 15.02.05	

function new_psi = psi_modif(psi,mod_p)

if psi.format ~= '01',
	disp('psi does not have the good format')
	return
end
if length(psi.psitbxfun.t) ~= 1,
	disp('the equilibrium must be found for one time only')
end

%isolate the divertor problem
a = plas_grid(psi);

% applies the transformation
switch mod_p(1)
case {0, 10}
	tmp = psi.psitbxfun.x;
	tmp(a==0 & tmp < 1) = 2-tmp(a==0 & tmp < 1);
	new_psi = psitbxpsi(tmp,psi.psitbxfun.grid,psi.psitbxfun.t,'01');
	new_psi.rmag = psi.rmag;
	new_psi.zmag = psi.zmag;
	new_psi.psimag = psi.psimag;
	
	
case {1,11}
	tmp = psi.psitbxfun.x;
	tmp(a==0 & tmp < 1) = 2-tmp(a==0 & tmp < 1);
	new_grid = psitbxgrid('Cylindrical','Grid',{psi.psitbxfun.grid.x{1}+mod_p(2),psi.psitbxfun.grid.x{2}});
	new_psi = psitbxpsi(tmp,new_grid,psi.psitbxfun.t,'01');
	new_psi.rmag = psi.rmag+mod_p(2);
	new_psi.zmag = psi.zmag;
	new_psi.psimag = psi.psimag;
	
case {2,12}
	%determines the points inside the plasma last closed shell
	tmp = psi.psitbxfun.x;
	tmp(a==0 & tmp < 1) = 2-tmp(a==0 & tmp < 1);
	[old_R,Z] = ndgrid(psi.psitbxfun.grid.x{1},psi.psitbxfun.grid.x{2});
	
	if mod_p(1) == 2 & mod_p(2) ~= 0,%check the possibility to use mod_p
	     tmp2 = psitbxp2p(psi,'FS');
	     if mod_p(2) < 0,
	     	if abs(mod_p(2)) > tmp2.psitbxfun.x(end,iround(tmp2.psitbxfun.grid.x{2},0))
			error('dR is too big in psi_modif')
		end
	     else 
	     	if mod_p(2) > tmp2.psitbxfun.x(end,end)
			error('dR is too big in psi_modif')
		end
	     end
	end
	new_R = old_R + mod_p(2)*(a.*(1-sqrt(tmp)));
	tmp = griddata(new_R,Z,tmp,old_R,Z,'cubic');
	if any(any(isnan(tmp))), disp('there are nans in 2'); end
	new_psi = psitbxpsi(tmp,psi.psitbxfun.grid,psi.psitbxfun.t,'01');
	new_psi.rmag = psi.rmag+mod_p(2);
	new_psi.zmag = psi.zmag;
	new_psi.psimag = psi.psimag;
	
case {3,13}
	tmp = psi.psitbxfun.x;
	tmp(a==0 & tmp < 1) = 2-tmp(a==0 & tmp < 1);
	new_grid = psitbxgrid('Cylindrical','Grid',{psi.psitbxfun.grid.x{1}+mod_p(2),psi.psitbxfun.grid.x{2}});
	tmp = psitbxpsi(tmp,new_grid,psi.psitbxfun.t,'01');
	tmp.rmag = psi.rmag+mod_p(2);
	tmp.zmag = psi.zmag;
	tmp.psimag = psi.psimag;
	psi = tmp;
	[old_R,Z] = ndgrid(psi.psitbxfun.grid.x{1},psi.psitbxfun.grid.x{2});
	
	if mod_p(1)== 3 & mod_p(3)~=0,
	      tmp = psitbxp2p(psi,'FS');
	      if mod_p(3) < 0,
	     	if abs(mod_p(3)) > tmp.psitbxfun.x(end,iround(tmp.psitbxfun.grid.x{2},0))
			error('dR is too big in psi_modif')
		end
	      else 
	     	if mod_p(3) > tmp.psitbxfun.x(end,end)
			error('dR is too big in psi_modif')
		end
	      end
	end
	
	new_R = old_R + mod_p(3)*(a.*(1-sqrt(psi.psitbxfun.x)));
	tmp = griddata(new_R,Z,psi.psitbxfun.x,old_R,Z,'cubic');
	if any(any(isnan(tmp))), disp('there are nans in 3'); end
	new_psi = psitbxpsi(tmp,psi.psitbxfun.grid,psi.psitbxfun.t,'01');
	new_psi.rmag = psi.rmag+mod_p(3);
	new_psi.zmag = psi.zmag;
	new_psi.psimag = psi.psimag;
		
case {4,14}
	tmp = psi.psitbxfun.x;
	tmp(a==0 & tmp < 1) = 2-tmp(a==0 & tmp < 1);
	new_grid = psitbxgrid('Cylindrical','Grid',{psi.psitbxfun.grid.x{1}+mod_p(2),psi.psitbxfun.grid.x{2}+mod_p(3)});
	new_psi = psitbxpsi(tmp,new_grid,psi.psitbxfun.t,'01');
	new_psi.rmag = psi.rmag + mod_p(2);
	new_psi.zmag = psi.zmag + mod_p(3);
	new_psi.psimag = psi.psimag;
otherwise
	error('Wrong input for fsd_mod.')
end

%to avoid the error with the fsd operation
new_psi = psitbxp2p(new_psi,'+0');
new_psi = psitbxp2p(new_psi,'01');
	
