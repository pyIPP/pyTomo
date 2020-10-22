function [data] = get_mpx(shot,t1,t2);
global_p;

% Function [data] = get_mpx(shot,t1,t2)
%
% The function loads the MPX data for discharges > 20500. rho is calculated
% as the tangential point between the lines of sight and the closest flux
% surface.
%
% Input:    shot         shot number
%           [t1 t2]   time interval to load
%
% Output:   
%           data      structured array with data and time base
%                     Channels are ordered from HFS to LFS
% History: 
%  24/10/01 new configuration included
% 
%                


sxrpath='/mac/camenen/matlab/mpx/';


disp('loading MPX data from tree')
mdsopen(shot)

if shot<19923
% taken from Sushkov's SXRVIEW routine + loads of guesswork
    % (see /home/sushkov/sxrview on HAL.)
    % I am not sure about the exact timing, the channels, as well as
    % tree structure.
    freq=2e5; dt=5e-6; n0=0;
    if shot>18093 
      dr = mdsdata('\atlas::dt100_001:slow:dec_rate');
      n0 = mdsdata('\atlas::dt100_001:slow:n0');
      np = mdsdata('\atlas::dt100_001:slow:n_points'); 
      sfreq=freq/dr; dt=dt*dr;
    else
      error('Does not work for shot#<18093');
      np=round((t2-t1)/dt);
      sfreq=freq;
    end
    t0=-.04+n0*5e-6;

    for i=1:64
      switch 1
      case shot<17900
        str=sprintf('\\atlas::dt100_001:channel_%.3d',i);
      case shot<17941
        str=sprintf('\\atlas::dt100_001:fast:channel_%.3d',i);
      case shot<18094
        str=sprintf('\\atlas::dt100_001:fast:channel_%.3d:raw',i);
      otherwise
        str=sprintf('\\atlas::dt100_001:slow:channel_%.3d:raw',i);
      end
      % sort channels
      if shot<18791
        if i<33
          n(i)=i*2-1;
        else
          n(i)=66-(i-32)*2;
        end
      else
        if i<33
          n(i)=66-i*2;
        else
          n(i)=(i-32)*2-1;
        end
      end
      data.data(:,n(i))=mdsdata(str);
		  end 
    [nt,nc]=size(data.data);
    data.dim  = {[t0:dt:t0+(nt-1)*dt]',[1:64]};
    % select data in defined time interval
    kk = find(data.dim{1}>=t1 & data.dim{1}<=t2);
    data.data=data.data(kk,:);
    data.dim{1}=data.dim{1}(kk);

%------------------------From Manini's file--------------------------------------------
elseif shot>19923 
% ordering from HFS to LFS   
  for ii=1:9
     eval(['data.data(:,65-(2*ii-1)) = mdsdata(''\ATLAS::DT100_NORTHEAST_001:CHANNEL_00' int2str(ii) ''');'])
  end

  for ii=10:32
     eval(['data.data(:,65-(2*ii-1)) = mdsdata(''\ATLAS::DT100_NORTHEAST_001:CHANNEL_0' int2str(ii) ''');'])
  end
  for ii=1:23
     eval(['data.data(:,65-2*ii) = mdsdata(''\ATLAS::DT100_NORTHEAST_002:CHANNEL_0' ...
     int2str(33-ii) ''');'])
  end
  for ii=24:32
     eval(['data.data(:,65-2*ii) = mdsdata(''\ATLAS::DT100_NORTHEAST_002:CHANNEL_00' ...
       int2str(33-ii) ''');'])
  end
  data.dim = {mdsdata('dim_of(\ATLAS::DT100_NORTHEAST_001:CHANNEL_001)'),[1:64]};
% time base incorrect, need to rescale
  if shot<20939
    data.dim{1}=(data.dim{1}+0.32)*0.9697-0.32;
  end;
% select data in defined time interval
  kk = find(data.dim{1}>=t1 & data.dim{1}<=t2);
  data.data=data.data(kk,:);
  data.dim{1}=data.dim{1}(kk);
else
  error(['Data not available for shot#',num2str(shot)]);
end
[nt,nc]=size(data.data);
% calibration 
load(sprintf('%smpx_calibration',sxrpath),'kk')
load(sprintf('%smpx_offset',sxrpath),'offset')
for i=1:64
  data.data(:,i)=data.data(:,i)*kk(i)-offset(i)*3.052e-4;
end
data.units = 'V';
data.dimunits ={'s', 'ch'};
mdsclose
% select data in defined time interval
kk = find(data.dim{1}>=t1 & data.dim{1}<=t2);
data.data=data.data(kk,:);
data.dim{1}=data.dim{1}(kk);
