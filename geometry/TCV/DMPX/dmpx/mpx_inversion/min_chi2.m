%min_chi2.m is used for the inversion of the DMPX stationary profile.
%	It finds the 1-D emissivity function g(rho) in flux coordinate.
%	It solves f = T*g with a least square fit. g has two edge 
%	conditions: dg/d(rho) = 0 at rho = 0 and g(1) = 0.
%
%
%[g,chi2,rho] = min_chi2(f,df,T,rho)
%
%
%Inputs:
%
%	f		data column vector (Nchords x 1)
%
%	df		absolute error on f (Nchords x 1)
%
%	T		transfer matrix (Nchords x Npixels)
%
%	rho		radial flux coordinates of the center of the pixels
%			(1 x Npixels). Warning: Npixels < Nchords !!!!
%
%
%Outputs:
%
%	g		emissivity profile
%
%	chi2		= sum((T*g - f).^2) / sum(df.^2); if chi2 > 1, the
%			result is completely false
%
%	rho		rho could eventually be modified.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Jonathan Rossel, 15.02.05
	
function [g,chi2,rho] = min_chi2(f,df,T,rho)


if length(rho) >= length(f),
	disp('there are too many pixels')
	return
end

%removes negative emissivities
f(f<0) = 0;

%initialisation
rho = sort(rho);
if rho(end) < 0.95,
	rho(end) = 0.95;
end

%constraints addition
C1 = [-1 1 zeros(1,length(rho)-2)]; %derivative close to rho = 0
C2 = [zeros(1,length(rho)-1), 1]; %value of g close to rho = 1
new_T = [T; C1; C2];
new_f = [f; 0; 0];
new_df = [df; 0; 0];

%weighting
tmp=(df<=0);
W(~tmp) = 1./df(~tmp);
W(tmp) = max(W(~tmp));
W=W';
W = diag([W; 10*max(W); 10*max(W)]);
new_T = W*new_T;
new_f = W*new_f;

%solving
g = new_T \ new_f;
chi2 = sum((T*g - f).^2) / sum(df.^2);
