%function mpxplot
%clear global
ww = which('mpxplot');
addpath([ww(1:end-9),'mpx_plot/']);
global_p;
shot=input('Enter shot number: ');
t_1=input('Enter t1: ');
t_2=input('Enter t2: ');
%calib=input('Calibration by default (0) or force the use of the new one (1) :');
%calib=0;

%mdsopen(shot);

detec='top';
if shot>=26609 
 I=input('Do you want top (1) or bottom (2) detector ?');
 if I==2, detec='bot'; end;
end 

if shot>=20030
	[mpx]=mpxdata(shot,'sg','time',[t_1 t_2],'detec',detec);
	eval(['signal=mpx.' detec '.signal;']);
else
	[signal]=get_mpx(shot,t_1,t_2);
	disp(['Warning shot<20030, you are using an old MPX routine that I did not check personaly' ])
end

% init array of windows names
windows_name(1,:)='window 1';
windows_name(2,:)='window 2';

% set current window number
wnd=1;
sws_3=0; zoo_m=0;
hold_w(1)=0; hold_w(2)=0;
n_hold(1)=1; n_hold(2)=1;

contr_val1=[1 0 0 0];
contr_val2=[1 0 0 0];

n_r=30; n_t=1; n_t_last=0;

xlab_2=3;

az_step=10; el_step=10;

% from load_mpx_data
 %%%%%%%%%%%%%%%%%%%% load psi %%%%%%%%%%%%%%%%%
if isempty(psi_data)==1
 [psi_data,psi_time,psi_R,psi_Z]=psidata(shot);
 psi_tmp=psi_data;
 %psi_time(length(psi_time))
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
 psi_1=psi_data(:,:,TT)';
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 load vessel
%if shot>=23807 
% load new_mpx_geometry
%else
% load mpx_geometry
%end
end

[n_tmax,n_rmax]=size(signal.data);

if shot>=20030
[n_tmax,n_rmax]=size(mpx.signal.data);
	eval(['chordZ(1,:)=mpx.' detec '.geom.wires.Z*100;']);
	eval(['chordR(1,:)=mpx.' detec '.geom.wires.R*100;']);
	chordZ(2,:)=repmat(76,1,n_rmax);
	chordR(2,:)=repmat(0,1,n_rmax);
	for ii=1:n_rmax,
	 chordR(:,ii)=interpos(11,[chordZ(1,ii) mpx.top.geom.slit.Z*100],[chordR(1,ii) mpx.top.geom.slit.R*100],chordZ(:,ii));
	end
else
	load mpx_geometry
end

t=signal.dim{1};


for n_tmp=1:n_rmax;
  r(n_tmp)=n_tmp;
end;  

n_rmin=1; n_tmin=1;
n_tmin1=1; n_tmax1=n_tmax;
n_rmin1=1; n_rmax1=n_rmax;

if exist('t')~=1; t=1:n_tmax; end;

if exist('r')~=1;
	if n_rmax>1;r=1:n_rmax;else; r(1)=1; end;
end;

if exist('viz')~=1;	
	for tmp_r=1:10;
	for tmp_t=1:10;
		vizZ(tmp_r)=1;
		vizR(tmp_r)=1;
		viz(tmp_t,tmp_r)=1;
	end;
	end;
end;
if exist('hordR')~=1;
	for r_tmp=1:n_rmax;
		hordR(1,r_tmp)=1; hordR(2,r_tmp)=2; hordZ(1,r_tmp)=1;
	        hordZ(2,r_tmp)=2;
	end;
end;

t_min=min(t); t_max=max(t);
r_min=min(r); r_max=max(r);
tmp=min(y);y_min=min(tmp);tmp=max(y);y_max=max(tmp);

n_t=round(length(t)/2);
n_r=round(length(r)/2);


init12w
y=signal.data;
plotter2(0)