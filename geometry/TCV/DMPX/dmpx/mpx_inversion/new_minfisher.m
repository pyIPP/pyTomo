%new_minfisher.m is used for the inversion of the DMPX stationary profile.
%	It finds the 1-D emissivity function g(rho) in flux coordinate.
%	It solves f = T*g with a minimization of the Fisher information. 
%	g has two edge conditions: dg/d(rho) = 0 at rho = 0 and g(1) = 0.
%
%
%[g,chi2,rho] = new_minfisher(f,df,T,rho)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	
	
function [g,chi2,rho] = new_minfisher(f,df,T,rho)

%removes negative emissivities
f(f<0) = 0;

%construction of the 1st order differential operator
rho = sort(rho);
if rho(end) < 0.95,
	rho(end) = 0.95;
end
delta_rho = diff(rho);
nrho = length(rho);
diff1 = zeros(nrho,nrho);
diff1(1,[1 2]) = [-1, 1] ./ delta_rho(1);
diff1(end,[nrho-1, nrho]) = [-1, 1] ./ delta_rho(end);
for ii = 2:nrho-1,
	diff1(ii,[ii,ii+1])=[-1,1]./2./delta_rho(ii); %best one to avoid oscillations in the profile
end


%normalisation
tmp=(df<=0);
norm_fac(~tmp)=1./df(~tmp);
tmp1=mean(norm_fac(~tmp));
norm_fac(tmp)=tmp1;
norm_fac=norm_fac';
T=T.*repmat(norm_fac,1,length(rho));
f=f.*norm_fac;
diff1=diff1./mean(df);


%weighting
W=diag([ones(size(df))]);


%solving
Wfish = eye(nrho,nrho);
H1 = diff1'*diff1;
T2 = T'*T;
niter = 0;
alpha_0=10*norm(T2,'fro')/norm(H1,'fro');
chi2target=0.1/length(f);

while 1, %loop for Wfish
	niter = niter+1;
	niter_a = 0;
	chi2_old=1e20;
	alpha = alpha_0;
	while 1, %loop for alpha
		niter_a = niter_a + 1;
		M = T2 + alpha*(H1);
		if rcond(M) < 2*eps,
			disp('Matrix badly conditionned')
			g = nan*zeros(size(T,2),1);
			chi2 = realmax;
			rho=NaN*zeros(1,size(T,2));
			return
		end
		
		g = M \ (T'*f);
		chi2 = mean((T*g - f).^2);
		
		%test
		if  niter_a > 100| ((chi2_old-chi2)/chi2_old<0.01 & niter_a>5),
			break
		else
			alpha=alpha/(1+abs(chi2-chi2target));
			chi2_old=chi2;
		end
	end
	old_Wfish = diag(Wfish);
	Wfish = 1 ./ g;
	tmp1=(g<0);
	tmp2=(g>0&g<=0.001*max(g));
	Wmax=max(Wfish); %maximum Wfish allowed
	Wfish(tmp1) = Wmax; 
	Wfish(tmp2) = mean(Wfish(~tmp1&~tmp2));
	Wfish(Wfish>Wmax)=Wmax;
	Wfish = diag(Wfish);
	Wfish = Wfish / norm(Wfish);
	
	test = mean((old_Wfish-diag(Wfish)).^2) / mean(old_Wfish.^2);
	if test < 0.01 | niter > 10,
		break
	end
	H1 = diff1'*Wfish*diff1;
end	
	
