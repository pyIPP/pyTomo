function [top_eff,bot_eff]=mpxeff(gasnum)

 % MPXEFF Calculates and plot both MPX detectors (top and bottom) efficiencies for the different
% absorber-holder positions
%
% SYNTAX
%
% [top_eff,bot_eff]=mpxeff(gas)
%
% INPUTS
%
% gas		= vector of detection gas number
% You can give one or several detection gas number in the vector:
% 1		= 90% Krypton + 10% Methan (CH4) detection gas (default)
% 2		= 90% Argon + 10% Methan (CH4) detection gas
% 3		= 90% Xenon + 10% Methan (CH4) detection gas
%
% OUTPUTS
%
% top_eff	= structure containing the TOP detector efficiency for
%		the 4 absorber-holder positions (i.e. for the 4 different
%		additional absorbers)
% bot_eff	= same structure for the BOTTOM detector
%
% EXAMPLES
%
% [top_eff,bot_eff]=mpxeff; % Only for Kr + CH4
% [top_eff,bot_eff]=mpxeff([1 2 3]); % For all gases
%
% Written by Loic Curchod (2008).
%

addpath /home/matlab/crpptbx-7.6.0/dmpx/mpx_efficiency

% Default detection gas

if ~exist('gasnum')
	gasnum = 1;
end

gas={'KrCHabs','ArCHabs','XeCHabs'};

for ii=gasnum

	%
	% TOP DETECTOR
	%

	figure;

	% Position 4 - No additional absorber

	detec = ['He229600, Be100 $ ' gas{ii} '8000'];
	pos4 = XraySignals(detec);
	pos4.detec=detec;
	semilogx(pos4.ev/1000,pos4.response);
	hold on;

	% Position 3 - 125 microns Be

	detec = ['He229600, Be225 $ ' gas{ii} '8000'];
	pos3 = XraySignals(detec);
	pos3.detec=detec;
	semilogx(pos3.ev/1000,pos3.response,'k');

	% Position 2 - 550 microns Be

	detec = ['He229600, Be650 $ ' gas{ii} '8000'];
	pos2 = XraySignals(detec);
	pos2.detec=detec;
	semilogx(pos2.ev/1000,pos2.response,'g');

	% Position 1 - 308 microns Al

	detec = ['He229600, Be100, Al308 $ ' gas{ii} '8000'];
	pos1 = XraySignals(detec);
	pos1.detec=detec;
	semilogx(pos1.ev/1000,pos1.response,'r');

	title(['DMPX TOP detector efficiency, absorb. ' gas{ii}(1:4) '4']);
	xlabel('E_{photon} [keV]');
	ylabel('absorption probability');
	legend('no additional abs.','125 \mum Be','550 \mum Be','308 \mum Al');
	grid off;
	axis([pos1.ev(1)/1000 pos1.ev(end)/1000 0 1]);
	set(gca,'XTicklabel',[1 10 100]);
	
	eval(['top_eff.' gas{ii} '.pos1=pos1;']);
	eval(['top_eff.' gas{ii} '.pos2=pos2;']);
	eval(['top_eff.' gas{ii} '.pos3=pos3;']);
	eval(['top_eff.' gas{ii} '.pos4=pos4;']);
	
	%
	% BOTH DETECTORS (first part: TOP detector)
	%
	
	hfboth=figure;
	
	semilogx(pos4.ev/1000,pos4.response);
	hold on;
	semilogx(pos3.ev/1000,pos3.response,'k');
	semilogx(pos2.ev/1000,pos2.response,'g');
	semilogx(pos1.ev/1000,pos1.response,'r');

	
	%
	% BOTTOM DETECTOR
	%

	figure;

	% Position 4 - No additional absorber

	detec = ['He229600, Be200, KrCH8000, AIR10000 $ ' gas{ii} '7600'];
	pos4 = XraySignals(detec);
	pos4.detec=detec;
	semilogx(pos4.ev/1000,pos4.response);
	hold on;

	% Position 3 - 125 microns Be
	
	detec = ['He229600, Be325, KrCH8000, AIR10000 $' gas{ii} '7600'];
	pos3 = XraySignals(detec);
	pos3.detec=detec;
	semilogx(pos3.ev/1000,pos3.response,'k');

	% Position 2 - 550 microns Be
	
	detec = ['He229600, Be750, KrCH8000, AIR10000 $ ' gas{ii} '7600'];
	pos2 = XraySignals(detec);
	pos2.detec=detec;
	semilogx(pos2.ev/1000,pos2.response,'g');
	
	% Position 1 - 308 microns Al
	
	detec = ['He229600, Be200, KrCH8000, AIR10000, Al308 $ ' gas{ii} '7600'];
	pos1 = XraySignals(detec);
	pos1.detec=detec;
	semilogx(pos1.ev/1000,pos1.response,'r');
	
	title(['DMPX BOTTOM detector efficiency, absorb. ' gas{ii}(1:4) '4']);
	xlabel('E_{photon} [keV]');
	ylabel('absorption probability');
	legend('no additional abs.','125 \mum Be','550 \mum Be','308 \mum Al');
	grid off;
	ax=axis;
	axis([pos1.ev(1)/1000 pos1.ev(end)/1000 0 1]);
	set(gca,'XTicklabel',[1 10 100]);

	eval(['bot_eff.' gas{ii} '.pos1=pos1;']);
	eval(['bot_eff.' gas{ii} '.pos2=pos2;']);
	eval(['bot_eff.' gas{ii} '.pos3=pos3;']);
	eval(['bot_eff.' gas{ii} '.pos4=pos4;']);
	
	%
	% BOTH DETECTORS (second part: BOTTOM detector)
	%
	
	figure(hfboth);
	hold on;
	semilogx(pos4.ev/1000,pos4.response,'--');
	semilogx(pos3.ev/1000,pos3.response,'--k');
	semilogx(pos2.ev/1000,pos2.response,'--g');
	semilogx(pos1.ev/1000,pos1.response,'--r');

	title(['DMPX detectors efficiency, absorb. ' gas{ii}(1:4) '4']);
	xlabel('E_{photon} [keV]');
	ylabel('absorption probability');
	legend('no add. abs. TOP','125 \mum Be TOP',...
	'550 \mum Be TOP','308 \mum Al TOP',...
	'no add. abs. BOT','125 \mum Be BOT',...
	'550 \mum Be BOT','308 \mum Al BOT');
	grid off;
	ax=axis;
	axis([pos1.ev(1)/1000 pos1.ev(end)/1000 0 1]);
	set(gca,'XTicklabel',[1 10 100]);
	
	%
	% BOTH DETECTORS (plotted in an other order)
	%
	
	figure;
	
	eval(['semilogx(top_eff.' gas{ii} '.pos4.ev/1000,top_eff.' gas{ii} '.pos4.response);']);
	hold on;
	eval(['semilogx(bot_eff.' gas{ii} '.pos4.ev/1000,bot_eff.' gas{ii} '.pos4.response,''--'');']);
	eval(['semilogx(top_eff.' gas{ii} '.pos3.ev/1000,top_eff.' gas{ii} '.pos3.response,''k'');']);
	eval(['semilogx(bot_eff.' gas{ii} '.pos3.ev/1000,bot_eff.' gas{ii} '.pos3.response,''--k'');']);
	eval(['semilogx(top_eff.' gas{ii} '.pos2.ev/1000,top_eff.' gas{ii} '.pos2.response,''g'');']);
	eval(['semilogx(bot_eff.' gas{ii} '.pos2.ev/1000,bot_eff.' gas{ii} '.pos2.response,''--g'');']);
	eval(['semilogx(top_eff.' gas{ii} '.pos1.ev/1000,top_eff.' gas{ii} '.pos1.response,''r'');']);
	eval(['semilogx(bot_eff.' gas{ii} '.pos1.ev/1000,bot_eff.' gas{ii} '.pos1.response,''--r'');']);

	title(['DMPX detectors efficiency, absorb. ' gas{ii}(1:4) '4']);
	xlabel('E_{photon} [keV]');
	ylabel('absorption probability');
	legend('no add. abs. TOP','no add. abs. BOT',...
	'125 \mum Be TOP','125 \mum Be BOT',...
	'550 \mum Be TOP','550 \mum Be BOT',...
	'308 \mum Al TOP','308 \mum Al BOT');
	grid off;
	ax=axis;
	eval(['axis([bot_eff.' gas{ii} '.pos1.ev(1)/1000 bot_eff.' gas{ii} '.pos1.ev(end)/1000 0 1]);']);
	set(gca,'XTicklabel',[1 10 100]);
end
