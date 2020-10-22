function [mpx_st,mpx_out] = mpxsawtooth(shot,t1,t2,std_fact_crash,std_fact_peak,dt,opt,plot_opt,mpx)

% function [mpx_st,mpx_out] = mpxsawtooth(shot,t1,t2,std_fact_crash,std_fact_peak,dt,opt,plot_opt,mpx)
%
% MPXSAWTOOTH analyzes the sawtooth activity based on DMPX data.
%
% INPUTS:
% shot		= shot number
% [t1,t2]	= time window for the analysis
% std_fact_crash = standard deviation factor for crashes detection.
%             Crashes with d(S)/dt < -1*std(d(S)/dt)*std_fact_crash are
%             considered as sawtooth crashes (where S is the DMPX signal).
%             Default: std_fact_crash = 0.5.
% std_fact_peak	= standard deviation factor for peaks detection.
%             Peaks with d(S)/dt > std(d(S)/dt)*std_fact_peak are
%             considered as peaks and withdrawn from the sawtooth crashes
%             list. Default = std_fact_peak = 1.75. Not used if negative.
% dt        = Takes profiles at t_crash +/- dt to find the inversion radius.
%             If empty or dt = 0, takes default times, i.e. times when the
%             signal first derivative crosses the detection threshold
%             around the crash times (t_crash_start and t_crash stop).
%             If length(dt) = 2, profiles are averaged on dt(2) before
%             looking for rho_inv. Time intervals for the profiles averaging
%             are indicated by color patches in the first figure. Examples:
%             dt = [dt1 dt2]: rho_inv is calculated using average profiles
%             in [t_crash-dt1-dt2:t_crash-dt1] and [t_crash+dt1:t_crash+dt1+dt2].
%             dt = [0 dt2]: rho_inv is calculated using average profiles in
%             [t_crash_start-dt2:t_crash_start] and [t_crash_stop:t_crash_stop+dt2]
% opt       = 1 or 11 or 21: loads only central DMPX chord and calculate sawtooth frequency only.
%             2 or 12 or 22: loads full DMPX signal to calculate both sawtooth frequency and inversion radius
%             If > 10, loads signal only in analysis time window [t1,t2] instead of [0,2] s.
%             If > 20, you need to provide DMPX inverted data in the mpx input structure.
%             Default: opt = 11.
% plot_opt	= 0: No plot.
%             1: Plots signal and analysis results.
%             2: Like 1, plus plots profiles around all crashes (in one
%             figure)
%             3: Like 2, but profiles in separate figures. (Careful, this
%             may create a large number of figures = number of sawtooth
%             crashes.)
%             Default: plot_opt = 1.
% mpx		= mpx structure from mpxdata.m containing data for the top detector.
%             Optional. If not given as input, mpx is loaded using mpxdata.m
%
% OUTPUTS:
% mpx_st	= output structure with fields:
% t_crash 	= sawtooth crashes times
% t_start_crash = sawtooth crashes start times
% t_stop_crash	= sawtooth crashes stop times
% tau_st	= sawtooth period
% f_st		= sawtooth frequency
% chan_inv	= estimated DMPX channel of sawtooth inversion
% rhopsi_inv	= corresponding rho_psi (of the flux surface which is tangent to the line of sight)
% rhovol_inv	= corresponding rho_vol (of the flux surface which is tangent to the line of sight)
% mean_prof_before_crash	= average profile on 0.2 ms before t_start_crash
% mean_prof_after_crash		= average profile on 0.2 ms after t_stop_crash
%
% DEFAULT CALL for sawtooth frequency analysis
% [mpx_st,mpx_out] = mpxsawtooth(shot,t1,t2,[],-1,[],1,1);
%
% EXAMPLES
% [mpx_st,mpx_out] = mpxsawtooth(31541,0.1,0.7,0.75,-1,[],12);	% ELM-free H-mode (with EBH at the edge)
% [mpx_st,mpx_out] = mpxsawtooth(31564,0.1,1.7,[],-1,[],12);	% An old standard shot
% [mpx_st,mpx_out] = mpxsawtooth(33026,0.7,1.9,1.75,-1,[],12);	% F. Felici's evil shot
% [mpx_st,mpx_out] = mpxsawtooth(38694,0.2,1.8,[],-1,[],12);	% A. Bortolon's ECH saturated ST
% [mpx_st,mpx_out] = mpxsawtooth(39635,0.1,1.6,[],-1,[],12);	% A recent standard shot
% [mpx_st,mpx_out] = mpxsawtooth(39670,0.2,1.85,0.7,-1,2.5e-3,12);	% A. Zhuchkova's ECH saturated ST
% [mpx_st,mpx_out] = mpxsawtooth(39669,1.4,1.6,0.7,-1,2.5e-3,12,2);	% A. Zhuchkova's ECH saturated ST

% Created by L. Curchod, CRPP, EPFL. 25.03.2010

if isempty(dt)
   dt_crash = 0;
   dt_mean = 0;
elseif length(dt)==1
    dt_crash = dt;
    dt_mean = 0;
else
    dt_crash = dt(1);
    dt_mean = dt(2);
end

if opt < 20
    fast_opt = 0;	% 0: 20 kHz data. 1: 200 kHz data.
    rhovol_opt = 1;	% 0: loads rhopsi only. 1: loads rhopsi and rhovol.
    detec = 'top';	% 'top' for top detector. 'bot' for bottom detector.
elseif opt > 20 & isfield(mpx.top.rho,'rhovol')
    rhovol_opt = 1;
    detec = 'top';
else
    rhovol_opt = 0;
    detec = 'top';
end

if detec == 'top'
    if opt < 20
        nchan = 64;	% Number of channels
        chan = 32;	% Channel number for sawtooth crash time analysis
    else
        nchan = size(mpx.top.signal.data,2);
        chan = 1;
    end
elseif detec == 'bot'
	nchan = 32;	% Number of channels
	chan = 16;	% Channel number for sawtooth crash time analysis
end

if ~exist('std_fact_crash','var');
	std_fact_crash = 0.5;
elseif isempty(std_fact_crash)
	std_fact_crash = 0.5;
end

if ~exist('std_fact_peak','var');
	std_fact_peak = 1.5;
elseif isempty(std_fact_peak)
	std_fact_peak = 1.5;
end

if ~exist('opt','var')
    opt = 11;
elseif isempty(opt)
    opt = 11;
elseif ~ismember(opt,[1 11 21 2 12 22])
    error('The ''opt'' input should be either 1, 11, 21, 2, 12 or 22.');  
end

if ~exist('plot_opt','var');
	plot_opt = 1;
elseif isempty(plot_opt);
	plot_opt = 1;
elseif ~ismember(plot_opt,[0 1 2 3])
	error('The ''plot_opt'' input should be either 0, 1, 2 or 3.');
end

if ~exist('mpx','var')
	%load_data = input('Are you sure you want to load DMPX data? [y/n] ','s');
	%if load_data=='n'
	%	return;
	%end
	if opt>10
        ta = t1;
        tb = t2;
    else
        if t1<0
            ta = t1;
        else
            ta = 0;
        end
        if t2>2
            tb = t2;
        else
            tb = 2;
        end
    end
    if opt == 1 | opt == 11
        [mpx] = mpxdata(shot,'svgf','time',[ta tb],'freq',fast_opt,'chords',chan,'vol',0);
    elseif opt == 2 | opt == 12
        [mpx] = mpxdata(shot,'svgfr','time',[ta tb],'freq',fast_opt,'detec',detec,'vol',rhovol_opt);
    end
	detec_struct = getfield(mpx,detec);
else
	if ~isfield(mpx,detec);
		error(['No ',detec,' field in the input mpx structure.']);
	else
		detec_struct = getfield(mpx,detec);
		if ~isfield(detec_struct,'signal')
			error(['No signal field in the input mpx.',detec,' structure.']);
		elseif size(detec_struct.signal.data,2)~=nchan & opt < 20
			error(['The signal.data matrix must have ',num2str(nchan),' channels.']);
		end
	end
end

mpx_out = mpx;

if 0
ok_chords = zeros(1,length(detec_struct.signal.dim{2}));
ok_chords(chan) = 1;
alpha = 10;
[signal_wo_spike] = remove_spike(detec_struct.signal,ok_chords,alpha);
end

iit = find(detec_struct.signal.dim{1}>=t1 & detec_struct.signal.dim{1}<=t2);
time = detec_struct.signal.dim{1}(iit);
data = detec_struct.signal.data(iit,chan);
profs = detec_struct.signal.data(iit,:);

%
% Plot calibrated data
%
if plot_opt
	hf = figure;
	pos = get(gcf,'Position');
	new_pos = pos;
	new_pos(4) = pos(4)+300;
	new_pos(2) = pos(2)-300;
	set(hf,'Position',new_pos);

	hs(1) = subplot(3,1,1); hold on; grid on;
	plot(detec_struct.signal.dim{1},detec_struct.signal.data(:,chan));
	plot(time,data,'k');
	xlabel('time [s]'); ylabel('Calibrated');
	title(['DMPX #',num2str(chan),' - Shot #',num2str(shot)]);
	ylim = get(gca,'YLim');
end

f_acq = 1/mean(diff(time));	% Acquisition frequency
N = 2;
if f_acq>100e3	% Bandpass between 20 Hz and 1000 Hz
	Wn = [0.0002 0.01];	% Bandpass between 1/10000 and 1/100 of f_acq
else
	Wn = [0.002 0.1];	% Bandpass between 1/1000 and 1/10 of f_acq
end
time = time(~isnan(data));
data = data(~isnan(data));
[B,A] = butter(N,Wn);
dataf = filtfilt(B,A,data);
%
% Plot filtered data
%
if plot_opt
	hs(2) = subplot(3,1,2); hold on; grid on;
	plot(time,dataf,'g');
	xlabel('time [s]'); ylabel('Filtered');
end

[time_prime,ii,jj] = unique(time);
[dataf_int,dataf_prime]=interpos(13,time_prime,dataf(ii),time_prime);
%
% Plot time derivative of filtered data
%
if plot_opt
	figure(hf);
	hs(3) = subplot(3,1,3); hold on; grid on;
	plot(time_prime,dataf_prime,'m');
	xlabel('time [s]'); ylabel('d(S_{filt})/dt');
end

max_tau_st = min(0.05,t2-t1); % 50 ms > longest sawtooth period ~ 20 ms
ntimes = floor(f_acq*max_tau_st);
nloops = max(floor(length(dataf_prime)/ntimes),1);
for ii=1:nloops
	if ii<nloops
		s(ii) = std(dataf_prime(1+(ii-1)*ntimes:ii*ntimes));	% Standard deviation in the interval
		s2(1+(ii-1)*ntimes:ii*ntimes,1) = ones(1,ntimes)*s(ii);	% Vector of st. dev.
	else
		s(ii) = std(dataf_prime(1+(ii-1)*ntimes:end));	% Standard deviation in the interval
		s2(1+(ii-1)*ntimes:length(dataf_prime),1) = ones(1,length(dataf_prime)-(nloops-1)*ntimes)*s(ii);
	end
end
%
% Plot standard deviation and ST detection threshold
%
if plot_opt
	plot(time_prime,s2,'--k');
	plot(time_prime,-s2,'--k');
	plot(time_prime,-s2*std_fact_crash,'-r');
    if std_fact_peak>0
        plot(time_prime,s2*std_fact_peak,'-r');
    end
	plot([time_prime(1),time_prime(end)],[0,0],'-.k');
end
%
% Find crashes times
%
iis = find(dataf_prime<= -s2*std_fact_crash);
diis = diff(iis);
%ii_start_crash = [iis(1);iis(find(diis>1)+1)]; % Don't take first crash and last crash
%ii_stop_crash = [iis(find(diis>1));iis(end)];
ii_start_crash = [iis(1);iis(find(diis>1)+1)];
ii_start_crash = ii_start_crash(2:end); % For safety, don't take first crash
ii_stop_crash = [iis(find(diis>1));iis(end)];
%
% Start with a start and stop with a stop
%
if ii_stop_crash(1)<=ii_start_crash(1)
	ii_stop_crash = ii_stop_crash(2:end);
end
if ii_start_crash(end)>=ii_stop_crash(end)
	ii_start_crash = ii_start_crash(1:end-1);
end
t_start_crash = time_prime(ii_start_crash);
t_stop_crash = time_prime(ii_stop_crash);
%
% Plot crashes starts and stops
%
if plot_opt
	plot(t_start_crash,dataf_prime(ii_start_crash),'or');
	plot(t_stop_crash,dataf_prime(ii_stop_crash),'sr');
end
iit_crash = zeros(length(ii_start_crash),1);
dataf_prime_crash = zeros(length(ii_start_crash),1);
for ii=1:length(ii_start_crash)
	ii_crash = [ii_start_crash(ii):ii_stop_crash(ii)];
	[dataf_prime_crash(ii),ii_min] = min(dataf_prime(ii_crash));
	iit_crash(ii) = ii_crash(ii_min);
end
t_crash = time_prime(iit_crash);
%
% Plot crashes times
%
if plot_opt
	plot(t_crash,dataf_prime_crash,'vr');
end
%
% Crashes period and frequency
%
tau_st = diff(t_crash);
f_st = 1./diff(t_crash);
%
% Find peaks times
%
if std_fact_peak>0
    iis = find(dataf_prime >= s2*std_fact_peak)
    diis = diff(iis);
    ii_start_peak = [iis(1);iis(find(diis>1)+1)];
    ii_stop_peak = [iis(find(diis>1));iis(end)];
    %
    % Start with a start and stop with a stop. Test
    %
    if 0
    if ii_stop_peak(1)<=ii_start_peak(1)
        ii_stop_peak = ii_stop_peak(2:end);
    end
    if ii_start_peak(end)>=ii_stop_peak(end)
        ii_start_peak = ii_start_peak(1:end-1);
    end
    end
    t_start_peak = time_prime(ii_start_peak);
    t_stop_peak = time_prime(ii_stop_peak);
    %
    % Plot peaks starts and stops
    %
    if plot_opt
        plot(t_start_peak,dataf_prime(ii_start_peak),'or');
        plot(t_stop_peak,dataf_prime(ii_stop_peak),'sr');
    end
    iit_peak = zeros(length(ii_start_peak),1);
    dataf_prime_peak = zeros(length(ii_start_peak),1);
    for ii=1:length(ii_start_peak)
        ii_peak = [ii_start_peak(ii):ii_stop_peak(ii)];
        [dataf_prime_peak(ii),ii_max] = max(dataf_prime(ii_peak));
        iit_peak(ii) = ii_peak(ii_max);
    end
    t_peak = time_prime(iit_peak);
    %
    % Plot peaks times
    %
    if plot_opt
        plot(t_peak,dataf_prime_peak,'vr');
    end
    t_crash_wo_peak = t_crash;
    iip = iround(t_crash,t_peak);
    t_crash_wo_peak(iip) = [];
    tau_st_wo_peak = diff(t_crash_wo_peak);
    f_st_wo_peak = 1./diff(t_crash_wo_peak);
end
%
% Plot tau_ST and f_ST
%
if plot_opt
	figure;
	pos = get(gcf,'Position');
	new_pos = pos;
	new_pos(4) = pos(4)+300;
	new_pos(2) = pos(2)-300;
	set(gcf,'Position',new_pos);
    %
    hs(4) = subplot(2,1,1); hold on; grid on;
    plot((t_crash(1:end-1)+t_crash(2:end))/2,tau_st*1e3);
    if std_fact_peak>0
        plot((t_crash_wo_peak(1:end-1)+t_crash_wo_peak(2:end))/2,tau_st_wo_peak*1e3,'--m');
	end
    xlabel('time [s]'); ylabel('\tau_{st} [ms]');
	title(['Shot #',num2str(shot)]);
	%
    hs(5) = subplot(2,1,2); hold on; grid on;
	plot((t_crash(1:end-1)+t_crash(2:end))/2,f_st);
    if std_fact_peak>0
        plot((t_crash_wo_peak(1:end-1)+t_crash_wo_peak(2:end))/2,f_st_wo_peak,'--m');
    end
	xlabel('time [s]'); ylabel('f_{st} [Hz]');
    figure(hf);
	subplot(3,1,1);
	for ii=1:length(t_crash)
		plot([t_crash(ii),t_crash(ii)],ylim,'--r');
    end
    if dt_crash==0 & dt_mean~=0
        for ii=1:length(t_crash)
            hp(1) = patch([t_start_crash(ii)-dt_mean t_start_crash(ii) t_start_crash(ii) t_start_crash(ii)-dt_mean],...
                [ylim(1) ylim(1) ylim(2) ylim(2)],'r');
            hp(2) = patch([t_stop_crash(ii) t_stop_crash(ii)+dt_mean t_stop_crash(ii)+dt_mean t_stop_crash(ii)],...
                [ylim(1) ylim(1) ylim(2) ylim(2)],'r');
            set(hp,'FaceAlpha',0.5,'LineStyle','-.','EdgeColor','r');
        end
    elseif dt_crash==0 & dt_mean==0
        for ii=1:length(t_crash)
            plot([t_start_crash(ii),t_start_crash(ii)],ylim,'-.r');
            plot([t_stop_crash(ii),t_stop_crash(ii)],ylim,'-.r');
        end 
    elseif dt_crash~=0 & dt_mean==0
       for ii=1:length(t_crash)
           plot([t_start_crash(ii),t_start_crash(ii)],ylim,'-.r');
           plot([t_stop_crash(ii),t_stop_crash(ii)],ylim,'-.r');
           plot([t_crash(ii)-dt_crash,t_crash(ii)-dt_crash],ylim,'-.m');
           plot([t_crash(ii)+dt_crash,t_crash(ii)+dt_crash],ylim,'-.m');
       end
    elseif dt_crash~=0 & dt_mean~=0
        for ii=1:length(t_crash)
            hp(1) = patch([t_crash(ii)-dt_crash-dt_mean t_crash(ii)-dt_crash t_crash(ii)-dt_crash t_crash(ii)-dt_crash-dt_mean],...
                [ylim(1) ylim(1) ylim(2) ylim(2)],'m');
            hp(2) = patch([t_crash(ii)+dt_crash t_crash(ii)+dt_crash+dt_mean t_crash(ii)+dt_crash+dt_mean t_crash(ii)+dt_crash],...
                [ylim(1) ylim(1) ylim(2) ylim(2)],'m'); 
            set(hp,'FaceAlpha',0.5,'LineStyle','-.','EdgeColor','m');
            plot([t_start_crash(ii),t_start_crash(ii)],ylim,'-.r');
            plot([t_stop_crash(ii),t_stop_crash(ii)],ylim,'-.r');
        end
    end
end
%
% INVERSION RADIUS ANALYSIS
%
if opt == 2 | opt == 12 | opt == 22
    
    if opt < 20
        chan_inv = zeros(length(t_crash),2);
        rhopsi_inv = zeros(length(t_crash),2);
        rhovol_inv = zeros(length(t_crash),2);
    else
        chan_inv = zeros(length(t_crash),1);
        rhopsi_inv = zeros(length(t_crash),1);
        rhovol_inv = zeros(length(t_crash),1);
    end
    mean_prof_before_crash = zeros(length(t_crash),nchan);
    mean_prof_after_crash = zeros(length(t_crash),nchan);
    if plot_opt==2
        hf = figure;
    end
    for ii=1:length(t_crash)
        if dt_crash==0
            t_start = t_start_crash(ii);
            t_stop = t_stop_crash(ii);
        else
            t_start = t_crash(ii)-dt_crash;
            t_stop = t_crash(ii)+dt_crash;
        end
        if dt_mean==0
            iit_before = iround(time,t_start);
            iit_after = iround(time,t_stop);
            mean_prof_before_crash(ii,:) = profs(iit_before,:);
            mean_prof_after_crash(ii,:) = profs(iit_after,:);
        else
            iit_before = find(time>t_start-dt_mean & time<t_start);
            iit_after = find(time>t_stop & time<t_stop+dt_mean);
            mean_prof_before_crash(ii,:) = mean(profs(iit_before,:),1);
            mean_prof_after_crash(ii,:) = mean(profs(iit_after,:),1);
        end
        % Test where profile before crash is above profile after crash
        test_chan = (mean_prof_before_crash(ii,:)>=mean_prof_after_crash(ii,:));
        % Find the profiles intersections
        dd_test_chan = diff(test_chan);
        if opt < 20
            ii_test_chan_LFS = find(dd_test_chan(32:end)==-1);
            ii_test_chan_HFS = find(dd_test_chan(1:32)==1);
        else
            ii_test_chan_LFS = find(dd_test_chan==-1);
            ii_test_chan_HFS = [];
        end
        if isempty(ii_test_chan_LFS)
            if opt < 20
                chan_inv(ii,2) = NaN;
            else
                chan_inv(ii) = NaN;
            end
            if isfield(detec_struct,'rho')
                iit_rho = iround(detec_struct.rho.time,t_crash(ii));
                if isfield(detec_struct.rho,'rhopsi');
                    rhopsi_inv(ii,2) = NaN;
                end
                if isfield(detec_struct.rho,'rhovol');
                    rhovol_inv(ii,2) = NaN;
                end
            end
        else
            if opt < 20
                chan_inv(ii,2) = 31+ii_test_chan_LFS(1);
                if isfield(detec_struct,'rho')
                    iit_rho = iround(detec_struct.rho.time,t_crash(ii));
                    if isfield(detec_struct.rho,'rhopsi');
                        rhopsi_inv(ii,2) = detec_struct.rho.rhopsi(iit_rho,chan_inv(ii,2));
                    end
                    if isfield(detec_struct.rho,'rhovol');
                        rhovol_inv(ii,2) = detec_struct.rho.rhovol(iit_rho,chan_inv(ii,2));
                    end
                end
            else
                chan_inv(ii) = ii_test_chan_LFS(1);
                if isfield(detec_struct,'rho')
                    iit_rho = iround(detec_struct.rho.time,t_crash(ii));
                    if isfield(detec_struct.rho,'rhopsi');
                        rhopsi_inv(ii) = detec_struct.rho.rhopsi(iit_rho,chan_inv(ii,1));
                    end
                    if isfield(detec_struct.rho,'rhovol');
                        rhovol_inv(ii) = detec_struct.rho.rhovol(iit_rho,chan_inv(ii,1));
                    end
                end
            end
            
        end
        if isempty(ii_test_chan_HFS) & opt < 20
            chan_inv(ii,1) = NaN;
            if isfield(detec_struct,'rho')
                iit_rho = iround(detec_struct.rho.time,t_crash(ii));
                if isfield(detec_struct.rho,'rhopsi');
                    rhopsi_inv(ii,1) = NaN;
                end
                if isfield(detec_struct.rho,'rhovol');
                    rhovol_inv(ii,1) = NaN;
                end
            end
        elseif  opt < 20
            chan_inv(ii,1) = 1+ii_test_chan_HFS(end);
            if isfield(detec_struct,'rho')
                iit_rho = iround(detec_struct.rho.time,t_crash(ii));
                if isfield(detec_struct.rho,'rhopsi');
                    rhopsi_inv(ii,1) = detec_struct.rho.rhopsi(iit_rho,chan_inv(ii,1));
                end
                if isfield(detec_struct.rho,'rhovol');
                    rhovol_inv(ii,1) = detec_struct.rho.rhovol(iit_rho,chan_inv(ii,1));
                end
            end
        end
        if plot_opt > 1
            if plot_opt == 2
                figure(hf); hold off;
            elseif plot_opt == 3
                figure; hold on;
            end
            if isfield(detec_struct.rho,'rhopsi')
                plot(detec_struct.rho.rhopsi(iit_rho,:),mean_prof_before_crash(ii,:),'b'); hold on;
                plot(detec_struct.rho.rhopsi(iit_rho,:),mean_prof_after_crash(ii,:),'r');
                if ~isempty(chan_inv(ii,:)) & sum(chan_inv(ii,:))~=0 & ~isnan(sum(chan_inv(ii,:)))
                    plot(detec_struct.rho.rhopsi(iit_rho,chan_inv(ii,:)),mean_prof_before_crash(ii,chan_inv(ii,:)),'ob');
                    plot(detec_struct.rho.rhopsi(iit_rho,chan_inv(ii,:)),mean_prof_after_crash(ii,chan_inv(ii,:)),'sr');
                end
                xlabel('\rho_{\psi}');
            else
                plot(mean_prof_before_crash(ii,:),'b'); hold on;
                plot(mean_prof_after_crash(ii,:),'r');
                if ~isempty(chan_inv(ii,:)) & sum(chan_inv(ii,:))~=0 & ~isnan(sum(chan_inv(ii,:)))
                    plot(chan_inv(ii,:),mean_prof_before_crash(ii,chan_inv(ii,:)),'ob');
                    plot(chan_inv(ii,:),mean_prof_after_crash(ii,chan_inv(ii,:)),'sr');
                end
                xlabel('DMPX channels');
            end
            ylabel('DMPX signal [a.u.]');
            title(['Crash at t = ',num2str(t_crash(ii)),' s']);
            legend('Before crash','After crash');
        end
    end
    %
    % Plot inversion radii
    %
    if plot_opt
        n_plot = 1;
        if isfield(detec_struct.rho,'rhopsi');
            n_plot = n_plot+1;
        end
        if isfield(detec_struct.rho,'rhovol');
            n_plot = n_plot+1;
        end
        figure;
        pos = get(gcf,'Position');
        new_pos = pos;
        new_pos(4) = pos(4)+300;
        new_pos(2) = pos(2)-300;
        set(gcf,'Position',new_pos);

        plot_num = 1;
        hs(5+plot_num) = subplot(n_plot,1,plot_num); hold on; grid on;
        plot(t_crash,chan_inv(:,1),'b');
        if opt < 20
            plot(t_crash,chan_inv(:,2),'r');
        end
        if opt<20
            ax = axis;
            axis([ax(1) ax(2) 1 64]);
        end
        xlabel('time [s]'); ylabel('MPX chan.');
        title(['Sawtooth inversion radii - Shot #',num2str(shot)]);
        if opt<20
            legend('HFS','LFS');
        end
        plot_num = plot_num+1;
        if isfield(detec_struct.rho,'rhopsi');
            hs(5+plot_num) = subplot(n_plot,1,plot_num); hold on; grid on;
            plot(t_crash,rhopsi_inv(:,1),'b');
            if opt < 20
                plot(t_crash,rhopsi_inv(:,2),'r');
            end
            ax = axis;
            if opt<20
                axis([ax(1) ax(2) -1 1]);
            else
                axis([ax(1) ax(2) 0 1]);
            end
            xlabel('time [s]'); ylabel('\rho_{\psi}');
            if opt<20
                legend('HFS','LFS');
            end
            plot_num = plot_num+1;
        end
        if isfield(detec_struct.rho,'rhovol');
            hs(5+plot_num)=subplot(n_plot,1,plot_num); hold on; grid on;
            plot(t_crash,rhovol_inv(:,1),'b');
            if opt < 20
                plot(t_crash,rhovol_inv(:,2),'r');
            end
            ax = axis; axis([ax(1) ax(2) -1 1]);
            xlabel('time [s]'); ylabel('\rho_{vol}');
            if opt<20
                legend('HFS','LFS');
            end
        end
    end
end
if plot_opt
    linkaxes(hs,'x');
end
%
% OUTPUT structure
%
mpx_st.t_crash = t_crash;
mpx_st.t_start_crash = t_start_crash;
mpx_st.t_stop_crash = t_stop_crash;
mpx_st.tau_st = tau_st;
mpx_st.f_st = f_st;
if std_fact_peak>0
    mpx_st.t_peak = t_peak;
    mpx_st.t_start_peak = t_start_peak;
    mpx_st.t_stop_peak = t_stop_peak;
end
%mpx_st.s2 = s2;
if opt == 2 | opt == 12 | opt == 22
    if opt<20
        mpx_st.chan_inv = chan_inv;
    else % For inverted data, no channel number any more
        mpx_st.ind_inv = chan_inv;
    end
    mpx_st.rhopsi_inv = rhopsi_inv;
    if rhovol_opt
        mpx_st.rhovol_inv = rhovol_inv;
    end
    mpx_st.mean_prof_before_crash = mean_prof_before_crash;
    mpx_st.mean_prof_after_crash = mean_prof_after_crash;
end
