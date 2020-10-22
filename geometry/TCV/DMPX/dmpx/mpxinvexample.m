% Example to run mpxinv.m:

shot=34081;time=1.3;dt=0.02;chords=[1:64];inv_method=3;fast=1;

% Inversion will be performed on [time-dt, time+dt]
% fast = 1 to use fast data (200 kHz), fast = 0 to use slow data (20 kHz)
% look in mpxinv.m help for inv_method

mpx=mpxdata(shot,'sgv','time',[time-dt time+dt],'freq',fast,'detec','top');
signal=mpx.top.signal;
sig_err=0.3*ones(size(signal.data)); %between 0.2 and 0.4
geom=mpx.top.geom;

psi=psitbxtcv(shot,time,'01');

ok_chords=zeros(1,length(signal.dim{2}));
ok_chords(chords)=1;

% First step (not required if you know fsd_mod) -> Find equilibrium modification
% for signal average on several tau_MHD or on a time interval during which plasma
% is stationnary

signal.data=mean(signal.data);
signal.dim{1}=mean(signal.dim{1});
sig_err=mean(sig_err);
fsd_mod=3;

[fsd_mod,new_psi,rho,g,chi2,T,I_chords,rho_HFS,g_HFS,chi2_HFS,T_HFS,I_chords_HFS,rho_LFS,g_LFS,chi2_LFS,T_LFS,I_chords_LFS]=...
	mpxinv(signal,sig_err,geom,psi,ok_chords,fsd_mod,inv_method);

% Compare inverted signal from LHF and HFS channels. Both profiles should be very
% similar at the edge.

figure;
plot(rho_HFS,g_HFS);
hold on;
plot(rho_LFS,g_LFS,'r');
xlabel('rho_{\psi}');
ylabel('Inverted signal [a.u.]');
legend('HFS chords','LFS chords');
title(['Shot #' num2str(shot) ', time=[' num2str(time-dt) ';' num2str(time+dt) ']']);

% You obtained an fsd_mod like fsd_mod=[3 -0.0091 -0.0082];
% You can save the mean profile inversion in a structure:

s34081_mean_inv = struct('I_chords',I_chords,'I_chords_HFS',I_chords_HFS,...
'I_chords_LFS',I_chords_LFS,'T',T,'T_HFS',T_HFS,'T_LFS',T_LFS,'chi2',chi2,...
'chi2_HFS',chi2_HFS,'chi2_LFS',chi2_LFS,'chords',chords,'dt',dt,'fast',fast,...
'fsd_mod',fsd_mod,'g',g,'g_HFS',g_HFS,'g_LFS',g_LFS,'geom',geom,'inv_method',...
inv_method,'new_psi',new_psi,'ok_chords',ok_chords,'psi',psi,'rho',rho,...
'rho_HFS',rho_HFS,'rho_LFS',rho_LFS,'shot',shot,'sig_err',sig_err,...
'signal',signal,'time',time);
save s34081_t1.28-1.32_mean_inv.mat s34081_mean_inv;

return;

% Second step -> inversion for each time

signal=mpx.top.signal;
sig_err=0.3*ones(1,size(signal.data,2)); % between 0.2 and 0.4

[fsd_mod,new_psi,rho,g,chi2,T,I_chords,rho_HFS,g_HFS,chi2_HFS,T_HFS,I_chords_HFS,rho_LFS,g_LFS,chi2_LFS,T_LFS,I_chords_LFS]=...
	mpxinv(signal,sig_err,geom,psi,ok_chords,fsd_mod,inv_method);

t=mpx.top.signal.dim{1}; % time

return;

% Store inverted data in a mpxdata-type structure.

mpx_invert.top.signal.data	= g';
mpx_invert.top.signal.dim{1} = t;
mpx_invert.top.signal.dim{2} = I_chords;
mpx_invert.top.rho.rhopsi = rho;
mpx_invert.top.rho.time = mean(t);

% Save the inversion data in a .mat file
s34081_inv = struct('chi2',chi2,'fsd_mod',fsd_mod,'inv_method',inv_method,...
    'mpx_invert',mpx_invert,'new_psi',new_psi,'psi',psi,'shot',shot,'T',T);
save s34081_t1.28-1.32_inv.mat s34081_inv;