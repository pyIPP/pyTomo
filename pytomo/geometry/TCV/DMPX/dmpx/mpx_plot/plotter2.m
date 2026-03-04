function plotter2(plot3d_);
global_p;

if shot==-1
 subplot(f1_1);
 plot(1,1)
 text(1,1,'No data for this shot')
 return
end
funct_var=100;

if n_t~=n_t_last;
 n_t_last=n_t; 
 if xlab_2==1|xlab_2==2   
  %psi_rho_chord=psitbxtcv(shot,t1,'01');
  %rho=rho_chord(chordR./100,chordZ./100,psi_rho_chord); %rho_psi
  if shot>=20030
  	mpx=mpxdata(shot,'r','rtime',t1);
  	eval(['rho=mpx.' detec '.rho.rhopsi;']);
  else %on triche un peu pour ne pas avoir a reecrire la routine mpxdata...pas bien :-( mais ca va plus vite :-)
  	mpx=mpxdata(20030,'g','detec','top'); 
    dcd=mpx.top.geom.dcd;
    psi=psitbxtcv(shot,t1,'01');
    rho=rhochords(dcd,psi);
  end
  if xlab_2==2,	
	%%%%%calcul de R pour Andrea
	fsd=psitbxtcv(shot,t1,'FS');
	R=rho*NaN;
	%HFS
	i_rho=find(rho<0),	
	tmp_rho=fliplr(abs(rho(i_rho))); %pour avoir les rho croissants, sinon psitbxgrid ne fonctionne pas
	fluxgrid=psitbxgrid('F','G',{tmp_rho,0,NaN},fsd);
    cylgrid_HFS=psitbxg2g(fluxgrid,'C');
 	R(i_rho)=flipud(cylgrid_HFS.x{1});
	%LFS
	i_rho=find(rho>=0);
  	tmp_rho=rho(i_rho);
	fluxgrid=psitbxgrid('F','G',{tmp_rho,-pi,NaN},fsd);
    cylgrid_LFS=psitbxg2g(fluxgrid,'C');
 	R(i_rho)=cylgrid_LFS.x{1};
  end	 
 end
 if sws_3==0 
  vis=1;
 end 
end

if xlab_2==1
 xlabel_2='\rho_{\psi}';
 r=rho;
end
if xlab_2==2
 xlabel_2='R [cm]';
 r=R';
end
if xlab_2==3
 xlabel_2='Channel';
 r=n_rmin:n_rmax;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(f1_1);
 
old_XLim1=get(gca,'XLim');
old_YLim1=get(gca,'YLim');
old_Xmode1=get(gca,'XLimMode');
old_Ymode1=get(gca,'YLimMode');

c_plot1=plot(t,y(:,n_r));

%set(gca,'YLim',YLim_tmp);
set(gca,'XLim',old_XLim1);
set(f1_1,'YLimMode','auto');
YLim_tmp=get(gca,'YLim');
  XL(1,:)=t(n_t)';
  XL(2,:)=XL(1,:);
 for i=1:length(n_t);
  YL(:,i)=YLim_tmp';
 end 
 t_line=line(XL,YL);
 clear XL YL
set(t_line,'LineStyle','-.');
		[tmp,n_tmin1]=min(abs(t-old_XLim1(1)));
		[tmp,n_tmax1]=min(abs(t-old_XLim1(2)));
		if n_tmin1<1;
			 n_tmin1=1;
		end;
	grid on	
 title (['#' num2str(shot)])
 xlabel 'Time [s]' 
 ylabel 'Isxr a.u.'         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subplot(f1_2);

		old_Xmode2=get(gca,'XLimMode');
		old_Ymode2=get(gca,'YLimMode');
		old_XLim2=get(gca,'XLim');
		old_YLim2=get(gca,'YLim');
  
	c_plot2=plot(r,y(n_t,:));
 YLim_tmp=get(gca,'YLim');
 XL(1,:)=r(n_r);
 XL(2,:)=XL(1,:);
 for i=1:length(n_r);
  YL(:,i)=YLim_tmp';
 end
 r_line=line(XL,YL);
 clear XL YL
 set(r_line,'LineStyle','-.');
%	set(gca,'YLim',YLim_tmp);
if old_Xmode2(1:1)=='m';
set(gca,'XLim',old_XLim2);
set(f1_2,'YLimMode','auto');
	[tmp,n_rmin1]=min(abs(r-old_XLim2(1)));
	[tmp,n_rmax1]=min(abs(r-old_XLim2(2)));
	if n_rmin1<1; 
		n_rmin1=1;
	end;
end

if exist('vis')==1
 visio
end 
if sws_3==0 
 subplot(f1_3)
 if isempty(v_lin)==0
  set(v_lin,'visible','off')
 end 
 ch_R(1,:)=chordR(1,n_r);
 ch_R(2,:)=chordR(2,n_r);
 ch_Z(1,:)=chordZ(1,n_r);
 ch_Z(2,:)=chordZ(2,n_r);
 v_lin=line(ch_R,ch_Z);
 set(v_lin,'visible','on')
end  
if sws_3==1
 v_lin=[];
end 

subplot(f1_2);
xlabel (xlabel_2) 
ylabel 'Isxr a.u.'
grid on

set(h_x_1,'String',sprintf('%g',t(n_t(1))));
set(h_y_1,'String',sprintf('%g',r(n_r(1))));

return;
   
