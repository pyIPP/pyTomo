function [psinorm_tbx,mpxgeom]=mpxlinesofsight(shot,time,numt,numt_mark,numb,numb_mark,ns)
%
% MPXLINESOFSIGHT plots selected lines-of-sight of the top and bottom wire chambers of the DMPX, on the top
% of a contour plot of normalized psi for a given shot at a given time
%
% SYNTAX
%
% [psinorm_tbx,mpxgeom]=mpxlinesofsight(shot,time,numt,numt_mark,numb,numb_mark)
%
% INPUTS
%
% shot		shot number
% time		time at which the normalized psi contour plot will be drawn
% numt		vector of the TOP detector lines-of-sight you want to plot
% numt_mark	vector of the TOP detector lines-of-sight you want to highlight
% numb		vector of the BOTTOM detector lines-of-sight you want to plot
% numb_mark	vector of the BOTTOM detector lines-of-sight you want to highlight
% ns        number of flux-surface contours to plot. Default = 11.
%
% if numt (or numb) = 0, all top (respect. bottom) lines-of-sight will be drawn
% if numt_mark (or numb_mark) < 0, no line-of-sight will be highlighted
%
% OUTPUTS
%
% psinorm_tbx	normalized psi as a psitbxtcv object
% mpxgeom	structure containing the top and bottom detectors geometry
%
% EXAMPLE
%
% [psinorm_tbx,mpxgeom]=mpxlinesofsight(35409,1.75,0,[35 40 47],0,[16 20 24]);
% [psinorm_tbx,mpxgeom]=mpxlinesofsight(35409,1.75,0,[35 40 47],[16 20 24],[],21);
%
% L.Curchod, July 2008. Based on fragments of mpxplot.m by Y.Camenen.
% Revised February 2009.
% Revised April 2009. Saved on the svn server.

if ~exist('ns')
    ns = 11;
end
if isempty(ns)
    ns = 11;
end

%addpath /home/matlab/crpptbx-7.6.0/dmpx/mpx_plot
%load vessel;

color={'b' 'r' [0 0.5 0] 'c' 'm' 'k'};

if numt==0
	numtop=[1:64];
else
	numtop=numt;
end

if numb==0
	numbot=[1:32];
else
	numbot=numb;
end

Z_stop_chords = 0.75;   % To stop chords at a given Z.

disp('Loading normalized psi for this shot...');
%mdsopen(shot);
%psi=tdi('\results::psi');
%mdsclose;
psinorm_tbx=psitbxtcv(shot,time,'01');	% Normalized psi from PSITBX
psinorm=psinorm_tbx.psitbxfun;

disp('Loading the DMPX geometry for this shot...');
if shot>=20030
	mpxgeom=mpxdata(shot,'g');
	chordZ(1,:)=mpxgeom.top.geom.wires.Z;
	chordR(1,:)=mpxgeom.top.geom.wires.R;
	chordZ(2,:)=repmat(Z_stop_chords,1,64); % To stop the chords
	for ii=1:64
		chordR(2,ii)=interpos(11,[chordZ(1,ii) mpxgeom.top.geom.slit.Z],[chordR(1,ii) mpxgeom.top.geom.slit.R],Z_stop_chords);
	end
else
	disp('Very old shot: loading an old geometry');
	load mpx_geometry;
	mpxgeom=0;
end

psi_R=psinorm.grid.x{1};
psi_Z=psinorm.grid.x{2};
psi_1=psinorm.x';

%psi_max=max(max(psi_1));
%psi_min=min(min(psi_1));

%V=(psi_max-psi_min)/20;
%VV=psi_min:V:psi_max;

VV = linspace(0,1,ns);

hf=figure;
set(hf,'Position',[500 100 400 600]);
tcvview('avt'); % Plot axes, vessel and tiles
tcvview('color','t',[0.75 0.75 0.75]); % Tiles in light gray
hold on; grid off;
contour(psi_R,psi_Z,psi_1,VV);
%hold on
%patch([Rv_in;Rv_in(1);Rv_out(1);flipud(Rv_out)], ...
%[Zv_in;Zv_in(1);Zv_out(1);flipud(Zv_out)],[0.5 0.5 0.5]);
axis equal;
axis([0.55 1.25 -0.85 0.90]);
xlabel('\itR\rm [m]');
ylabel('\itZ\rm [m]');
hc=colorbar;
%hcxlabel=get(hc,'Xlabel');
%set(hcxlabel,'String','\psi_n [Tm^2]');
%set(hcxlabel,'Position',[0.3696 0 1]);
hcylabel=get(hc,'Ylabel');
set(hcylabel,'String','\psi_n');
for ii=setdiff(numtop,numt_mark)
	hl=line([chordR(1,ii) chordR(2,ii)],[chordZ(1,ii) chordZ(2,ii)],...
    'Color','b','LineWidth',0.25);
end
if numt_mark>=0
	for ii=1:length(numt_mark)
		hl=line([chordR(1,numt_mark(ii)) chordR(2,numt_mark(ii))],...
		[chordZ(1,numt_mark(ii)) chordZ(2,numt_mark(ii))],'Color',color{mod(ii-1,6)+1});
		set(hl,'LineWidth',2);
		ht=text(chordR(2,numt_mark(ii)),chordZ(2,numt_mark(ii)),num2str(numt_mark(ii)),'Color',color{mod(ii-1,6)+1});
		%ht=text(chordR(2,numt_mark(ii)),0.83,num2str(numt_mark(ii)),'Color',color{mod(ii-1,6)+1});
		set(ht,'BackgroundColor',[1 1 1],'FontSize',14);
	end
end
if numt==0
	titre='DMPX, top lines of sight';	
else
	titre=['DMPX, top lines of sight [' num2str(numtop) ']'];
end
title(titre);
ht=text(0.6,0.84,['#' num2str(shot) ', \itt\rm = ' num2str(time) ' s']);
set(ht,'FontSize',12);

%
% Plots lines of sight of bottom detector
%

if shot>26554 %isfield(mpxgeom,'bot')

	chordZbot(1,:)=mpxgeom.bot.geom.wires.Z;
	chordRbot(1,:)=mpxgeom.bot.geom.wires.R;
	chordZbot(2,:)=repmat(Z_stop_chords,1,32); % To stop the chords
	for ii=1:32
		chordRbot(2,ii)=interpos(11,[chordZbot(1,ii) mpxgeom.bot.geom.slit.Z],...
		[chordRbot(1,ii) mpxgeom.bot.geom.slit.R],Z_stop_chords);
	end

	hf=figure;
    set(hf,'Position',[500 100 400 600]);
	tcvview('avt'); % Plot axes, vessel and tiles
    tcvview('color','t',[0.75 0.75 0.75]); % Tiles in light gray
    hold on; grid off;
    contour(psi_R,psi_Z,psi_1,VV);
    %hold on;
    %patch([Rv_in;Rv_in(1);Rv_out(1);flipud(Rv_out)], ...
	%[Zv_in;Zv_in(1);Zv_out(1);flipud(Zv_out)],[0.5 0.5 0.5]);
	axis equal;
	axis([0.55 1.25 -0.85 0.90]);
	xlabel('\itR\rm [m]');
	ylabel('\itZ\rm [m]');
	hc=colorbar;
	%hcxlabel=get(hc,'Xlabel');
	%set(hcxlabel,'String','\psi_n [Tm^2]');
	%set(hcxlabel,'Position',[0.3696 -0.01 1]);
    hcylabel=get(hc,'Ylabel');
	set(hcylabel,'String','\psi_n');
	for ii=setdiff(numbot,numb_mark)
		hl=line([chordRbot(1,ii) chordRbot(2,ii)],[chordZbot(1,ii) chordZbot(2,ii)],'Color','m');
		set(hl,'LineWidth',0.25);
	end
	if numb_mark>=0
		for ii=1:length(numb_mark)
			hl=line([chordRbot(1,numb_mark(ii)) chordRbot(2,numb_mark(ii))],...
			[chordZbot(1,numb_mark(ii)) chordZbot(2,numb_mark(ii))],'Color',color{mod(ii-1,6)+1});
			set(hl,'LineWidth',2);
			ht=text(chordRbot(2,numb_mark(ii)),chordZbot(2,numb_mark(ii)),num2str(numb_mark(ii)),'Color',color{mod(ii-1,6)+1});
			set(ht,'BackgroundColor',[1 1 1],'FontSize',14);
		end
	end
	if numb==0
		titre='DMPX, bottom lines of sight';	
	else
		titre=['DMPX, bottom lines of sight [' num2str(numbot) ']'];
	end
	title(titre);
	ht=text(0.6,0.84,['#' num2str(shot) ', \itt\rm = ' num2str(time) ' s']);
	set(ht,'FontSize',12);

	hf=figure;
    set(hf,'Position',[500 100 400 600]);
	tcvview('avt'); % Plot axes, vessel and tiles
    tcvview('color','t',[0.75 0.75 0.75]); % Tiles in light gray
    hold on; grid off;
    contour(psi_R,psi_Z,psi_1,VV);
	%hold on;
	%patch([Rv_in;Rv_in(1);Rv_out(1);flipud(Rv_out)], ...
	%[Zv_in;Zv_in(1);Zv_out(1);flipud(Zv_out)],[0.5 0.5 0.5]);
	axis equal;
	axis([0.55 1.25 -0.85 0.90]);
	xlabel('\itR\rm [m]');
	ylabel('\itZ\rm [m]');
	hc=colorbar;
	%hcxlabel=get(hc,'Xlabel');
	%set(hcxlabel,'String','\psi_n [Tm^2]');
	%set(hcxlabel,'Position',[0.3696 -0.01 1]);
    hcylabel=get(hc,'Ylabel');
	set(hcylabel,'String','\psi_n');
	for ii=numtop
		hl=line([chordR(1,ii) chordR(2,ii)],[chordZ(1,ii) chordZ(2,ii)],'Color','b');
		set(hl,'LineWidth',0.25);
	end
	for ii=numbot
		hl=line([chordRbot(1,ii) chordRbot(2,ii)],[chordZbot(1,ii) chordZbot(2,ii)],'Color','m');
		set(hl,'LineWidth',0.25);
	end
	titre='DMPX lines of sight';
	title(titre);
	ht=text(0.60,0.84,['#' num2str(shot) ', \itt\rm = ' num2str(time) ' s']);
	set(ht,'FontSize',12);

	hf=figure;
    set(hf,'Position',[500 100 880 600]);
	subplot(1,2,1);
	tcvview('avt'); % Plot axes, vessel and tiles
    tcvview('color','t',[0.75 0.75 0.75]); % Tiles in light gray
    hold on; grid off;
    contour(psi_R,psi_Z,psi_1,VV);
	%hold on;
	%patch([Rv_in;Rv_in(1);Rv_out(1);flipud(Rv_out)], ...
	%[Zv_in;Zv_in(1);Zv_out(1);flipud(Zv_out)],[0.5 0.5 0.5]);
	axis equal;
	axis([0.55 1.25 -0.85 0.90]);
	xlabel('\itR\rm [m]');
	ylabel('\itZ\rm [m]');
    hc=colorbar;
	%hcxlabel=get(hc,'Xlabel');
	%set(hcxlabel,'String','\psi_n [Tm^2]');
	%set(hcxlabel,'Position',[0.3696 -0.01 1]);
    hcylabel=get(hc,'Ylabel');
	set(hcylabel,'String','\psi_n');
	for ii=setdiff(numtop,numt_mark)
		hl=line([chordR(1,ii) chordR(2,ii)],[chordZ(1,ii) chordZ(2,ii)],'Color','b');
		set(hl,'LineWidth',0.25);
	end
	if numt_mark>=0
		for ii=1:length(numt_mark)
			hl=line([chordR(1,numt_mark(ii)) chordR(2,numt_mark(ii))],...
			[chordZ(1,numt_mark(ii)) chordZ(2,numt_mark(ii))],'Color',color{mod(ii-1,6)+1});
			set(hl,'LineWidth',2);
			ht=text(chordR(2,numt_mark(ii)),chordZ(2,numt_mark(ii)),num2str(numt_mark(ii)),'Color',color{mod(ii-1,6)+1});
			set(ht,'BackgroundColor',[1 1 1],'FontSize',14);
		end
	end
	if numt==0
		titre='DMPX, top lines of sight';	
	else
		titre=['DMPX, top lines of sight [' num2str(numtop) ']'];
	end
	title(titre);
	ht=text(0.60,0.84,['#' num2str(shot) ', \itt\rm = ' num2str(time) ' s']);
	set(ht,'FontSize',12);

	subplot(1,2,2);
    tcvview('avt'); % Plot axes, vessel and tiles
    tcvview('color','t',[0.75 0.75 0.75]); % Tiles in light gray
    hold on; grid off;
    contour(psi_R,psi_Z,psi_1,VV);
    %hold on;
	%patch([Rv_in;Rv_in(1);Rv_out(1);flipud(Rv_out)], ...
	%[Zv_in;Zv_in(1);Zv_out(1);flipud(Zv_out)],[0.5 0.5 0.5]);
	axis equal;
	axis([0.55 1.25 -0.85 0.90]);
	xlabel('\itR\rm [m]');
	ylabel('\itZ\rm [m]');
	hc=colorbar;
	%hcxlabel=get(hc,'Xlabel');
	%set(hcxlabel,'String','\psi_n [Tm^2]');
	%set(hcxlabel,'Position',[0.3696 -0.01 1]);
    hcylabel=get(hc,'Ylabel');
	set(hcylabel,'String','\psi_n');
	for ii=setdiff(numbot,numb_mark)
		hl=line([chordRbot(1,ii) chordRbot(2,ii)],[chordZbot(1,ii) chordZbot(2,ii)],'Color','m');
		set(hl,'LineWidth',0.25);
	end
	if numb_mark>=0
		for ii=1:length(numb_mark)
			hl=line([chordRbot(1,numb_mark(ii)) chordRbot(2,numb_mark(ii))],...
			[chordZbot(1,numb_mark(ii)) chordZbot(2,numb_mark(ii))],'Color',color{mod(ii-1,6)+1});
			set(hl,'LineWidth',2);
			ht=text(chordRbot(2,numb_mark(ii)),chordZbot(2,numb_mark(ii)),num2str(numb_mark(ii)),'Color',color{mod(ii-1,6)+1});
			set(ht,'BackgroundColor',[1 1 1],'FontSize',14);
		end
	end
	if numb==0
		titre='DMPX, bottom lines of sight';	
	else
		titre=['DMPX, bottom lines of sight [' num2str(numbot) ']'];
	end
	title(titre);
	ht=text(0.60,0.84,['#' num2str(shot) ', \itt\rm = ' num2str(time) ' s']);
	set(ht,'FontSize',12);
end
