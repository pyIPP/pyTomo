function [psi,psi_time,psi_R,psi_Z]=psidata(shotin);
global_p;
shot = shotin;
%mdsopen('frank::',shot);
mdsopen(shot);
psi=mdsdata('\results::psi').*sign(mean(mdsdata('\results::i_p')));
psi_time=mdsdata('dim_of(\results::psi,2)');
psi_R=mdsdata('dim_of(\results::psi)')*100;
psi_Z=mdsdata('dim_of(\results::psi,1)')*100;
%%%%%psi_toff=psi(:,:,Toff);
mdsclose

