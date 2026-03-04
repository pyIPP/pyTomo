%Removes the spikes (arcs in the wire chamber) of the DMPX signal
%From a routine by Jonathan Rossel (diploma work 2004)
%
%	function [signal_wo_spike] = remove_spike(signal,ok_chords,alpha);
%Inputs:
%	signal:	dim :time*chords
%		ex:	mpx=mpxdata(28357,'s');
%			signal=mpx.top.signal;
%	ok_chords: array with 1 if the chord is ok, 0 if not (default: ok_chords=ones(1,length(signal.dim{2})); )
%	alpha:	parameter for spike detection (default=10)
%
% Changelog:
%
%	14.04.09	J.Rossel	Adapt griddata call to 7.x

function [signal_wo_spike] = remove_spike(signal,ok_chords,alpha);

%inputs check and initialisation
if exist('signal')~=1
 error('remove_spike:WrongInput',['Input ''signal'' required'])
end 
if isfield(signal,'dim')&isfield(signal,'data')
 data=signal.data;
 time=signal.dim{1};
else
 error('remove_spike:WrongInput',['Input ''signal'' has not the required fields'])
end
if exist('ok_chords')~=1, 
 ok_chords=ones(1,length(signal.dim{2}));
elseif isempty(find(length(ok_chords)==size(data)))
 error('remove_spike:WrongInput',['Input ''signal'' and ''ok_chords'' must have one dimension in common']);
end
if exist('alpha')~=1, alpha=10; end 

max_length_spike=5; %max spike width (nb of chords) before 2D interpolation

%check the presence of NaN in the signal
I=find(isnan(mean(data))); %index of chords with at least one NaN
test=any(~isnan(data(:,I))); %chords with NaNs and real values
if sum(test)>0,
  warning(['Lonely NaNs in the signal for chord(s): ' num2str(I(test)) '. Chord(s) removed'])     	
end
ok_chords(I)=0;
disp(['Chords removed: ' num2str(find(~ok_chords))])
I_ok=find(ok_chords);
N_chords=sum(ok_chords);
N_time=length(time);
signal_wo_spike=signal;
signal_wo_spike.data(:,~ok_chords)=NaN;

%version check required by griddata (14.04.09)
verstr = version;
majorver = str2num(verstr(1));
if majorver >= 7, 
	blnMatlab7 = true;
else
	blnMatlab7 = false;
end

%%%%%%%%%%%%%% spike detection %%%%%%%%%%%%%%%%%%%
disp('Spike detection')		       
%detect the top of spikes time derivative greater than mean_der+alpha*std_der
der_data=abs(diff([zeros(1,N_chords);data(:,I_ok)],1,1)./repmat(diff([0;time]),1,N_chords));
n_points=min(2000,length(time)); %average the data by groups of n_points (in the time direction)
I=n_points-rem(N_time,n_points); %number of lines to be added to der_data
der_data_bis=reshape([der_data;repmat(NaN,I,N_chords)],n_points,(N_time+I)/n_points,N_chords);
mean_der_data_bis=repmat(mean(der_data_bis),[n_points 1 1]);
mean_der_data=reshape(mean_der_data_bis,N_time+I,N_chords);
mean_der_data(end-I+1:end,:)=[];
std_der_data_bis=repmat(std(der_data_bis),[n_points 1 1]);
std_der_data=reshape(std_der_data_bis,N_time+I,N_chords);
std_der_data(end-I+1:end,:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%new version: both last blocs are calculated together (remove the problem of end spikes)
if I > 0,
	end_length = n_points + rem(N_time,n_points);
	tmp = der_data(end-end_length+1:end,:);
	mean_der_data(end-end_length+1:end,:) = repmat(mean(tmp),end_length,1);
	std_der_data(end-end_length+1:end,:) = repmat(std(tmp),end_length,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = find(der_data >(mean_der_data+alpha*std_der_data));
no_spike=true([N_time N_chords]);
no_spike(index)=false;

%extend the spike zone until the time derivative changes sign
diff2s_data=abs(diff([zeros(1,N_chords);sign(diff([zeros(1,N_chords);data(:,I_ok)],1,1))],1,1)); %change of sign of the diff along time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%New version> flat zone can happen at the top of a spike (saturation)...
%diff2s_data(find(diff2s_data == 1))=2; %if there is a flat zone, it is also considered as the end of the spike
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diff2s_data(index)=3; %contains 3 inside top of spikes
for ii = 1:N_chords,
	tmp = find(diff2s_data(:,ii) == 2); %index of change of sign or flat zone outside the spikes
	tmp2 = find(diff2s_data(:,ii) == 3); %index of top of the spikes
	tmp2 = unique(ifloor(tmp,tmp2)); %use unique in case there are consecutive 3's
	tmp3 = find(tmp2 == 0); % happen if no change of sign at the beginning of the first spike
	tmp4 = tmp2+1;
	tmp2(tmp3)=[];
	spf = [ones(1,length(tmp3)), tmp(tmp2)']; %first index of the spikes
	tmp3 = find(tmp4 > length(tmp)); % happen if no change of sign at the end of the last spike
	tmp4(tmp3)=[];
	spl = [tmp(tmp4)',N_time*ones(1,length(tmp3))]; %last index of the spikes
	tmp = [];
	for jj = 1:length(spf),
		tmp = [tmp,spf(jj):spl(jj)];
	end
	no_spike(tmp,ii)=false;
end 
%look for large spike in the chords direction
%attention: find is columnwise, thus '
diff_no_spike=diff([ones(N_time,1) no_spike],1,2)'; %chords*time, -1 spike begins, 1 spike ends
[I_sp_beg,J_sp_beg]=find(diff_no_spike==-1); 
I=find(sum(diff_no_spike,1)~=0); %columns where a spike begins but does not end
diff_no_spike(end,I)=1; %ends the spike
[I_sp_end,J_sp_end]=find(diff_no_spike==1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%new version to have an exact correspondence -1/first zero and 1/last zero
if ~isempty(I),
	tmp = any(repmat(J_sp_end,1,length(I)) == repmat(I,length(J_sp_end),1),2); %vector, 1 if J_sp_end is the index of a corrected column
	tmp = I_sp_end == N_chords & tmp;  %1 if I_sp_end is the chord index of a corrected spike end
else
	tmp = false(1,length(I_sp_end));
end
I_sp_end(~tmp) = I_sp_end(~tmp) - 1;
len_sp=I_sp_end-I_sp_beg + 1; %spike length in the chord direction
I_long_spike = len_sp>max_length_spike;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%len_sp=I_sp_end-I_sp_beg; %spike length in the chord direction
%I_long_spike= find(len_sp>max_length_spike);
%short_spike=~no_spike;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%new version adds
[I_sp_beg_short,I_sp_end_short,J_sp_beg_short] = ...
	deal(I_sp_beg(~I_long_spike),I_sp_end(~I_long_spike),J_sp_beg(~I_long_spike));
[I_sp_beg_long,I_sp_end_long,J_sp_beg_long] = ...
	deal(I_sp_beg(I_long_spike),I_sp_end(I_long_spike),J_sp_beg(I_long_spike));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for ii=1:length(I_long_spike),
%	short_spike(J_sp_beg(I_long_spike(ii)),I_sp_beg(I_long_spike(ii)):I_sp_end(I_long_spike(ii)))=false;
%end


%%%%%%%%%%%%%% spike removal %%%%%%%%%%%%%%%%%%%%%
disp('Spike interpolation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%new version> add a boundary interpolation. Useful for small deviation
first_ch = data(:,I_ok(1));
last_ch = data(:,I_ok(end));
%interpolates first and last chord along time, but will reinterpolated along chords
%needed because bound. cond. must be inside the input interval for interpos
tmp2 = 1:N_time;
if any(~no_spike(:,1)),
	if length(find(~no_spike(:,1))) == 1,
		tmp = [tmp2(~no_spike(:,1)),tmp2(~no_spike(:,1)),tmp2(~no_spike(:,1))];
		tmp = interpos(13,tmp2(no_spike(:,1)),first_ch(no_spike(:,1)),tmp,10);
		tmp = tmp(1);
	else
		tmp = interpos(13,tmp2(no_spike(:,1)),first_ch(no_spike(:,1)),tmp2(~no_spike(:,1)),10);
	end
	first_ch(~no_spike(:,1)) = tmp;
end
if any(~no_spike(:,end)),
	if length(find(~no_spike(:,end))) == 1,
		tmp = [tmp2(~no_spike(:,end)),tmp2(~no_spike(:,end)),tmp2(~no_spike(:,end))];
		tmp = interpos(13,tmp2(no_spike(:,end)),first_ch(no_spike(:,end)),tmp,10);
		tmp = tmp(1);
	else
		tmp = interpos(13,tmp2(no_spike(:,end)),first_ch(no_spike(:,end)),tmp2(~no_spike(:,end)),10);
	end
	last_ch(~no_spike(:,end)) = tmp;
end
bc_HFS = zeros(N_time,1); %used for boundary conditions
tmp = no_spike(:,1) & no_spike(:,2);
bc_HFS(tmp) = diff(data(tmp,I_ok([1,2])),1,2)./(I_ok(2)-I_ok(1));
if ~isempty(tmp2(~tmp)),
	if length(tmp2(~tmp)) == 1,
		tmp3 = [tmp2(~tmp), tmp2(~tmp), tmp2(~tmp)];
		tmp3 = interpos(13,find(tmp),bc_HFS(tmp),tmp3,20);
		tmp3 = tmp3(1);
	else
		tmp3 = interpos(13,find(tmp),bc_HFS(tmp),tmp2(~tmp),20);
	end
	bc_HFS(~tmp) = tmp3;
end
bc_LFS = zeros(N_chords,1); %used for boundary conditions
tmp = no_spike(:,end-1) & no_spike(:,end);
bc_LFS(tmp) = diff(data(tmp,I_ok([N_chords,N_chords-1])),1,2)./(I_ok(end)-I_ok(end-1));
if ~isempty(tmp2(~tmp)),
	if length(tmp2(~tmp)) == 1,
		tmp3 = [tmp2(~tmp), tmp2(~tmp), tmp2(~tmp)];
		tmp3 = interpos(13,find(tmp),bc_LFS(tmp),tmp3,20);
		tmp3 = tmp3(1);
	else
		tmp3 = interpos(13,find(tmp),bc_LFS(tmp),tmp2(~tmp),20);
	end
	bc_LFS(~tmp) = tmp3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%interpolation in the chord direction for small spikes
%I=find(any(short_spike==true,2)); %indices of lines where there is at least 1 short spike
I = unique(J_sp_beg_short); %indices of lines where there is at least 1 short spike
for ii=1:length(I), %for each time where there is a small spike in the profile
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%new version (old one has end problems)
	%removes one additional chord before and after arc (improves the interpolation) if no other spikes
	J_sp_ind = find(J_sp_beg_short == I(ii));
	ind_interp = []; %index of chords which will be interpolated
	if length(J_sp_ind) == 1 & I_sp_beg_short(J_sp_ind) > 2 & I_sp_end_short(J_sp_ind) < N_chords-1,
		ind_interp = I_ok(I_sp_beg_short(J_sp_ind)-1:I_sp_end_short(J_sp_ind)+1);
	else
		for jj = 1:length(J_sp_ind),
			ind_interp = [ind_interp,I_ok(I_sp_beg_short(J_sp_ind(jj)):I_sp_end_short(J_sp_ind(jj)))];
		end
	end
	ind_for_interp = intersect(setdiff(I_ok,ind_interp),[I_ok(1),I_ok(end),I_ok(no_spike(I(ii),:))]);
    if length(ind_for_interp)<3
        ind_for_interp = [1 2 3 64];
    elseif length(ind_for_interp)<4 %interpos needs at least 4 input points (cubic spline)
        if ind_for_interp(end-1)~=63
            ind_for_interp = [ind_for_interp(1:end-1),ind_for_interp(end-1)+1,ind_for_interp(end)];
        else
            ind_for_interp = [ind_for_interp(1:end-2),ind_for_interp(end-1)-1,ind_for_interp(end-1:end)];
        end
    end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	
	%ind_interp=I_ok(short_spike(I(ii),:)); %index of chords which will be interpolated
	%tmp=diff([1 no_spike(I(ii),:)]); %removes one additional chord before and after arc (improves the interpolation)
	%i_m1=find(tmp==-1); %debut de spike
	%i_1=find(tmp==1);	%fin de spike
	%no_sp=no_spike(I(ii),:);
	%if ~isempty(i_m1)&(I_ok(max(1,i_m1-1))==I_ok(i_m1)-1),  
	%	for jj=1:length(i_m1)
	%		if (i_m1(jj)-2>=1)&no_sp(i_m1(jj)-2)~=0, %do not remove chord if it is alone (tmp(im_1-1)==1)
	%			no_sp(i_m1(jj)-1)=0; 
	%		end;
	%	end;
	%end
	%if ~isempty(i_1)&(I_ok(i_1-1)==I_ok(i_1)-1), 
	%	for jj=1:length(i_1)-1
	%		if no_sp(i_1(jj)+1)~=0, %do not remove chord if it is alone 
	%			no_sp(i_1(jj))=0; 
	%		end;
	%	end;
	%end	
	%ind_for_interp=I_ok(no_sp); %index of chords used to interpolate
	
	if length(ind_interp)==1,
		ind_interp=[ind_interp ind_interp ind_interp]; %interpos must have at least two values for xout
		%clear tmp
		tmp=interpos(13,ind_for_interp,[first_ch(I(ii)),data(I(ii),ind_for_interp(2:end-1)),last_ch(I(ii))],...
				ind_interp,10,[1,1],[bc_HFS(I(ii)),bc_LFS(I(ii))]);
		%tmp=interpos(13,ind_for_interp,data(I(ii),ind_for_interp),ind_interp,10);
		signal_wo_spike.data(I(ii),ind_interp)=tmp(1);
	else
		%signal_wo_spike.data(I(ii),ind_interp)=interpos(13,ind_for_interp,
		%data(I(ii),ind_for_interp),ind_interp,10);
        tmp = interpos(13,ind_for_interp,...
				[first_ch(I(ii)),data(I(ii),ind_for_interp(2:end-1)),last_ch(I(ii))],...
				ind_interp,10,[1,1],[bc_HFS(I(ii)),bc_LFS(I(ii))]);
        signal_wo_spike.data(I(ii),ind_interp)=tmp;

	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%new version: no need to interpolate twice...
	no_spike(I(ii),ifloor(I_ok,ind_interp)) = true;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%interpolation in the chord (whole profile) and in time direction (+ and - N points) for big spikes
ind_t=unique(J_sp_beg_long); %time index for which there is a long spike
for ii=1:length(ind_t),	
	
	N=10;
	pb=true;
	while pb==true&N<38,
		ind_dt=max(ind_t(ii)-N,1):min(ind_t(ii)+N,N_time);
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%New version: find the middle chord of the large spikes and smooth it along
		%time, keeps its data. Then, griddata
		tmp = find(J_sp_beg_long == ind_t(ii));
		for jj = 1:length(tmp), %loop on the long spikes at ind_t(ii)
			chord_ind = floor(mean([I_sp_beg_long(tmp(jj)),I_sp_end_long(tmp(jj))]));
			ind_for_interp = find(no_spike(ind_dt,chord_ind));
			ind_interp = find(~no_spike(ind_dt,chord_ind));
			if length(ind_interp) > 1 & length(ind_for_interp) > 2,
				tmp5 = find(ind_dt(ind_interp) == ind_t(ii)); %ind of ind_interp that is in the good spike
				%look for consecutive index from tmp
				tmp2 = diff(ind_interp);
				%beginnings of consecutive sequences
				tmp3 = find(tmp2 > 1)+1; 
				if tmp2(1)==1, tmp3 = [1,tmp3]; end %add the first index
				%ends of consecutive sequences
				tmp4 = find(tmp2 > 1);
				if tmp2(end) == 1, tmp4 = [tmp4,length(ind_interp)]; end %add the last index
				%redefine ind_interp
				ind_interp = ind_interp(tmp3(max(find(tmp3 <= tmp5))):tmp4(min(find(tmp4 >= tmp5))));
				if length(ind_interp) < 3,
					continue
				else
					signal_wo_spike.data(ind_dt(ind_interp),I_ok(chord_ind)) = ...
						interpos(13,ind_for_interp,signal_wo_spike.data(ind_dt(ind_for_interp),I_ok(chord_ind)),...
							ind_interp,20);
					no_spike(ind_dt(ind_interp),chord_ind) = true;
				end
			else
				continue
			end
		end
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%new version: expression of X and Z
		X = repmat(ind_dt',1,N_chords);
		X = X(no_spike(ind_dt,:));
		Z = signal_wo_spike.data(ind_dt,I_ok);
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		%X=repmat(ind_dt,1,N_chords);
		%X=X(no_spike(ind_dt,:))';
		Y=repmat(I_ok,length(ind_dt),1);
		Y=Y(no_spike(ind_dt,:));
		%Z=data(ind_dt,I_ok);
		Z=Z(no_spike(ind_dt,:));
		XI=ind_t(ii);
		YI=I_ok(~no_spike(ind_t(ii),:))'; %new version> column vector
		%if blnMatlab7,
			%tmp=griddata(X,Y,Z,XI,YI,'cubic',{'Qt','Qbb','Qc','Qz'}); %14.04.09: addition to avoid problems with griddata (in 7.x).
																	%'Qt','Qbb','Qc' are std options. 'Qz' adds a point to reduce round off error
																	%hopefully, it just makes the calculation longer.
		%else
			tmp=griddata(X,Y,Z,XI,YI,'cubic');
		%end
		if isempty(find(isnan(tmp))), 
			pb=false;	
		end
		N=N+5;
	end
	if ~isempty(find(isnan(tmp))),
		disp(['Pb in the large spike interpolation at t=' num2str(time(ind_t(ii))) ' for chords ' num2str(YI')]); 
	else
		signal_wo_spike.data(XI,YI)=tmp;
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%new version: use the already interpolated values for the next ones
		no_spike(XI,ifloor(I_ok,YI)) = true;
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	end	
	clear tmp X Y Z XI YI ind_dt
end

if any(any(~no_spike)),
	error('Some detected spikes were not removed')
end

disp([num2str(length(len_sp)) ' spikes removed (more or less) succesfully'])
	
	
	
