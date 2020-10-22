%SndOrdReg.m is used for the inversion of the DMPX stationary profile.
%	It finds the 1-D emissivity function g(rho) in flux coordinate.
%	It solves f = T*g with a minimization of the second order 
%	derivative. g has two edges conditions: dg/d(rho) = 0 at 
%	rho = 0 and g(1) = 0.
%
%
%[g,chi2,rho] = SndOrdReg(f,df,T,rho,chi2target)
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
%			(1 x Npixels)
%
%
%Optional inputs:
%
%	chi2target	chi2 reached at the end of the inversion.
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
%See also:
%
%	min_chi2.m 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Jonathan Rossel, 15.02.05
	

function [g,chi2,rho] = SndOrdReg(f,df,T,rho,chi2target)

if nargin < 5 | isempty(chi2target),
	chi2target = [];
end

%removes negative emissivities
f(f<0) = 0;

%find the target chi2
nrho = length(rho);
if isempty(chi2target),
	tmp = floor(nrho / length(f));
	tmp =  1:2*tmp:nrho;
	tmp2 = diff(tmp);
	rho_modif = rho(tmp);
	T_modif = zeros(length(f),length(rho_modif));
	for ii = 1:length(rho_modif)-1,
		T_modif(:,ii) = sum(T(:,[tmp(ii):tmp(ii)+tmp2(ii)-1]),2); 
	end
	T_modif(:,end) = sum(T(:,tmp(end):end),2);
	[tmp,chi2,tmp] = min_chi2(f,df,T_modif,rho_modif);
	chi2target = min([(1-chi2)/2+chi2,2*chi2,1]);
	chi2target = max(chi2target,0.1);%the 0.1 term comes from the necessity to have smoothing
end

%construction of the 2nd order differential operator
rho = sort(rho);
if rho(end) < 0.95,
	rho(end) = 0.95;
end
tmp = diff(rho);
tmp = filter([0.5 0.5],1,tmp);
delta_rho2 = (tmp(2:end)).^2;
diff2 = zeros(length(rho)-2,length(rho));
for ii = 1:size(diff2,1),
	diff2(ii,ii) = 1/delta_rho2(ii);  
	diff2(ii,ii+1) = -2/delta_rho2(ii);
	diff2(ii,ii+2) = 1/delta_rho2(ii);
end

%constraints on g
C1 = [-1 1 zeros(1,length(rho)-2)]; %derivative close to rho = 0
C2 = [zeros(1,length(rho)-1), 1]; %value of g close to rho = 1
new_T = [T; C1; C2];
new_f = [f; 0; 0];
new_df = [df; 0; 0];

%weighting
W = 1./df;
tmp = find(isinf(W));
W(tmp) = max(W(~isinf(W)));
W = diag([W; 10*max(W); 10*max(W)]);
new_T = W*new_T;
new_f = W*new_f;

%solving
T2 = new_T'*new_T;
D2 = diff2'*diff2;
alpha = trace(T2)/trace(D2);
niter = 0;
while 1, %loop for alpha, the regularization parameter
	niter = niter + 1;
	M = T2 + alpha*D2;
	if rcond(M) < 2*eps,
		disp('the matrix is badly conditionned, SndOrdReg')
		g = nan*zeros(size(T,2),1);
		chi2 = realmax;
		return
	end
	g = M \ (new_T'*new_f);

	%chi2
	chi2 = sum((T*g - f).^2) / sum(df.^2);
	
	%test
	if abs(chi2-chi2target) <= 0.005 | niter > 50,
		break
	else
		alpha = alpha/(chi2+(1-chi2target));
	end
end

%test the smoothing of g
d_g = diff(g);
extr_g_ind = find(diff(sign(d_g)) == 2 | diff(sign(d_g)) == -2)+1; %local extrema
tmp = diff(extr_g_ind); %tmp = 1 if 2 extrema are consecutive;
tmp = find(tmp==1);
tmp = unique([tmp;tmp+1]);
extr_g_ind = extr_g_ind(tmp); %indices of consecutive extrema
if isempty(extr_g_ind),
	return
end
tmp = [];
for ii = 1:length(extr_g_ind),
	tmp(ii) = mean(abs(d_g([extr_g_ind(ii)-1,extr_g_ind(ii)])));
end
if any(tmp > 0.5*std(d_g)), %if the profile is smooth, the variations must be very small
	disp('there are large oscillations, chi2target is increased')
	if chi2target == 1,
		disp('the errors on the MPX signal are not sufficient')
		df = df*2;
		[g,chi2,rho] = SndOrdReg(f,df,T,rho,0.1);
	else
		chi2target = min(2*chi2target,1);
		[g,chi2,rho] = SndOrdReg(f,df,T,rho,chi2target);
	end
end



