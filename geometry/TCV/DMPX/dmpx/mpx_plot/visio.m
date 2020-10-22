function visio
global_p;

t1=t(n_t(1));
if psi_time(length(psi_time))>t1
 tmp=find(psi_time>t1);
 if tmp(1)>1
  TT=tmp(1)-1; clear tmp
 else
  TT=1;
 end 
else 
 TT=(length(psi_time));
end  
psi_1=psi_tmp(:,:,TT)';

%rtmp=find(max(psi_1)>=0);
%ztmp=find(max(psi_1')>=0);
%if rtmp(1)>1
% psi_Rmin=psi_R(rtmp(1)-1);
%else 
% psi_Rmin=psi_R(rtmp(1));
%end 

%psi_Rmax=psi_R(rtmp(length(rtmp))+1);
%if ztmp(1)>1
% psi_Zmin=psi_Z(ztmp(1)-1);
%else
% psi_Zmin=psi_Z(ztmp(1));
%end 
%psi_Zmax=psi_Z(ztmp(length(ztmp))+1);

psi_max=max(max(psi_1));
V=psi_max/15;
VV=0:V:psi_max;

subplot(f1_2);
tmp=get(gca,'XLim');
[mm,nn]=find(r>=tmp(1));
nn1=nn(1);
[mm,nn]=find(r<=tmp(2));
nn2=nn(length(nn));
 
subplot(f1_3);
contour(psi_R,psi_Z,psi_1,VV);
patch([Rv_in;Rv_in(1);Rv_out(1);flipud(Rv_out)], ...
[Zv_in;Zv_in(1);Zv_out(1);flipud(Zv_out)],[0.5 0.5 0.5])  

for r_tmp=nn1:nn2;
	v_lin3(r_tmp)=line([chordR(1,r_tmp) chordR(2,r_tmp)],[chordZ(1,r_tmp) chordZ(2,r_tmp)]);
	set(v_lin3(r_tmp),'LineStyle',':');
	com_tmp=sprintf('set(v_lin3(r_tmp),''ButtonDownFcn'','';n_r=%i;plotter2(0);'')',r_tmp);
	eval(com_tmp)
end;

axis equal 
axis([min(Rv_out) max(Rv_out) min(Zv_out) max(Zv_out)]);
grid on
drawnow
v_lin=[];

return
