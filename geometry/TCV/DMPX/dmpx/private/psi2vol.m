function rho_vol=psi2vol(shot,time,rho_psi)
global_p;

mdsopen(shot);
vol=tdi('\results::psitbx:vol');
I=iround(vol.dim{2},time);
Vmax=max(vol.data(:,I));
x=sqrt(vol.data(:,I)./Vmax);


xout=[0:0.0001:1];
[yout]=interpos(2,vol.dim{1}',x,xout);
I=iround(xout,rho_psi);
rho_vol=yout(I);
