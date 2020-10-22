function [axuvstruct]=axuv_get_axuv(shotno,tbegin,tend,varargin)
%[axuvstruct]=axuv_get_axuv(shotno,tbegin,tend,varargin)
%e. g.:
%	[axuvstruct]=axuv_get_axuv(45103);
%	[axuvstruct]=axuv_get_axuv(45104,0.2,0.3);
%	[axuvstruct]=axuv_get_axuv(45109,[],[],'shotno_offset',45103);


%varargin:
%	'shotno_offset':	shot number if the data acquisition started after 0

% modifs: 2013.03.19: workaround for LOS#31 (BL)
shotno_offset=axuv_varargvalue('shotno_offset',varargin);

display('AXUV data download has started.');

if ~exist('tbegin')
    tbegin=[];
end
if ~exist('tend')
    tend=[];
end

if isstr(tbegin) || isstr(tend)
    mdsdisconnect;
    error('tbegin and tend must be a number or an empty variable.');
end


bolomap(:,1)=[1 20 2 19 3 18 4 17 5 16 6 15 7 14 8 13 9 12 10 11];
bolomap(:,2)=[43 62 44 61 45 60 46 59 47 58 48 57 49 56 50 55 51 54 52 53];
bolomap(:,3)=[1 94 2 93 3 92 4 91 5 90 6 89 7 88 8 87 9 86 10 85];
bolomap(:,4)= [43 42 44 41 45 40 46 39 47 38 48 37 49 36 50 35 51 34 52 33];
bolomap(:,5)=[85 84 86 83 87 82 88 81 89 80 90 79 91 78 92 77 93 76 94 75];
bolomap(:,6)=[21 42 22 41 23 40 24 39 25 38 26 37 27 36 28 35 29 34 30 33];
bolomap(:,7)=[65 84 66 83 67 82 68 81 69 80 70 79 71 78 72 77 73 76 74 75];

boards(:,[1:2 4:7])=repmat([11 11 12 12 13 13],20,1);
boards(:,3)=[12 11 12 11 12 11 12 11 12 11 12 11 12 11 12 11 12 11 12 11];

if shotno>=44499  %the cabling changed after the recommissioning!!!
    bolomap(:,[3 4 5])=flipud(bolomap(:,[3 4 5]));
    boards(:,[3 4 5])=flipud(boards(:,[3 4 5]));
    axuvstruct.geometry=axuv_calc_los(axuv_detectorpos('axuv'));
    if exist('axuv_finiteetend_axuv.mat')
        load axuv_finiteetend_axuv
        axuvstruct.geometry.etend=etendstruct.etend;
        etendstruct=[];
    end
    gain=400e3;
elseif shotno<44499 && shotno>=34666
    axuvstruct.geometry=axuv_calc_los(axuv_detectorpos('axuv_old'));
    gain=20e3;
else
    error('The routine is not suitable for downloading these data.');
end

%mdsopen('frank::tcv_shot',shotno);
mdsopen(shotno);
t0=mdsvalue('\ATLAS::DT196_AXUV_001:STOP_TRIG');	%triggering time
per=mdsvalue('slope_of(\ATLAS::DT196_AXUV_001:EXT_CLOCK_IN)');	%acquisition period (maybe sampling time)
ns=mdsvalue('\ATLAS::DT196_AXUV_001:MAX_SAMPLES')*1000; 	%ns is number of samples acquired
mdsdisconnect;

time_full=(t0+(0:(ns-1))*per)';
if time_full(1)<0
    tminoffsetindex=find(time_full<0,1,'first');
    tmaxoffsetindex=find(time_full<0,1,'last');
elseif time_full(end)>=2
    tminoffsetindex=find(time_full>2,1,'first');
    tmaxoffsetindex=find(time_full>2,1,'last');
else
    tmaxoffsetindex=[]; tminoffsetindex=[];
end
% if isempty(tbegin) && isempty(tend)
%     time=time_full;
%     tbeginindex=[];
%     tendindex=[];
% else
if isempty(tbegin)
    tbegin=time_full(1);
elseif time_full(1)>tbegin
    warning(['Data acquisition started at ' num2str(time_full(1)) ' sec.']);
    tbegin=(time_full(1));
end
if isempty(tend)
    tend=time_full(end);
elseif time_full(end)<tend
    warning(['Data acquisition finished at ' num2str(time_full(end)) ' sec.']);
    t_end=time_full(end);
end
indexes=find(time_full>=tbegin & time_full<=tend);
tbeginindex=indexes(1);
tendindex=indexes(end-1);
time=time_full(indexes(2:end));
% end


data_raw=zeros(length(time),140);
mdsopen('crpppc250.epfl.ch::axuv,$',shotno);
if ~isempty(tmaxoffsetindex)
    for i=1:140
        [data_raw(:,i),offset(i)]=get_chord_offset(boards(i),bolomap(i),tbeginindex,tendindex,tminoffsetindex,tmaxoffsetindex,shotno);
    end
    if ~isempty(shotno_offset)
        warning('Offset was measured in this shot. shotno_offset value is ignored.')
        shotno_offset=shotno;
    end
    mdsdisconnect;
else
    if ~isempty(shotno_offset)
        for i=1:140
            [data_raw(:,i)]=get_chord(boards(i),bolomap(i),tbeginindex,tendindex);
        end
        mdsdisconnect;
        mdsopen('frank::tcv_shot',shotno_offset);
        t0_for_offset_other=mdsvalue('\ATLAS::DT196_AXUV_001:STOP_TRIG');	%triggering time
        per_for_offset_other=mdsvalue('slope_of(\ATLAS::DT196_AXUV_001:EXT_CLOCK_IN)');	%acquisition period (maybe sampling time)
        ns_for_offset_other=mdsvalue('\ATLAS::DT196_AXUV_001:MAX_SAMPLES')*1000; 	%ns is number of samples acquired
        mdsdisconnect;
        time_for_offset_other=(t0_for_offset_other+(0:(ns_for_offset_other-1))*per_for_offset_other)';
        tminoffsetindex_other=1;
        tmaxoffsetindex_other=max(find(time_for_offset_other<0));
        if isempty(tmaxoffsetindex_other)
            mdsdisconnect;
            error('The shot which should have been used to calculate the offset has no offset.');
        end
        mdsopen('crpppc250.epfl.ch::axuv,$',shotno_offset);
        for i=1:140
            [offset(i)]=get_offset(boards(i),bolomap(i),tminoffsetindex_other,tmaxoffsetindex_other);
        end
        mdsdisconnect;
    else
        mdsdisconnect;
        error(['No offset data because acquisition started at ' num2str(t0) ' sec. Use an offset from a different shot with the shotno_offset option.']);
    end
end

axuvstruct.data=zeros(size(data_raw));
for i=1:size(data_raw,2)
    axuvstruct.data(:,i)=(data_raw(:,i)-offset(i))/0.24/gain/axuvstruct.geometry.etend(i);
end
% 10/2^15:	[Volt/bit] conversion factor of the ADC
% 0.24:		[A/W] conversion factor according to Boivin's article

axuvstruct.time=time;
axuvstruct.raw.data=data_raw;
axuvstruct.raw.offset=offset;
axuvstruct.raw.gains=repmat(gain,1,140);
axuvstruct.raw.convfact=0.24;
axuvstruct.raw.shotno_offset=shotno_offset;
axuvstruct.name='axuv';
axuvstruct.shotno=shotno;
display('AXUV data download has finished.');

function [b,offset]=get_chord_offset(board, channel,tbeginindex,tendindex,tminoffsetindex,tmaxoffsetindex,shotno)
if isempty(tbeginindex)
    b=mdsdata(['.BOARD_0' int2str(board) ':CHANNEL_0' num2str(channel,'%0.2d')])*10/2^15;
    if (size(b,1) ~= 1 && size(b,1) ~= 0)
        offset=mean(b(tminoffsetindex:tmaxoffsetindex));
    end
else
    b=mdsdata(['.BOARD_0' int2str(board) ':CHANNEL_0' num2str(channel,'%0.2d') '[' int2str(tbeginindex) ':' int2str(tendindex) ']'])*10/2^15;
    if (size(b,1) ~= 1 && size(b,1) ~= 0)
        offset=mean(mdsdata(['.BOARD_0' int2str(board) ':CHANNEL_0' num2str(channel,'%0.2d')...
            '[' int2str(tminoffsetindex) ':' int2str(tmaxoffsetindex) ']']))*10/2^15;
    end
end
%    disp(['Board ',num2str(board),' Channel ',num2str(channel)])
if shotno>=45523&&board==11&&channel==48
    b(b<-5)=b(b<-5)+10;
    if isempty(tbeginindex)
        if (size(b,1) ~= 1 && size(b,1) ~= 0)
            offset=mean(b(tminoffsetindex:tmaxoffsetindex));
        end
    else
        if (size(b,1) ~= 1 && size(b,1) ~= 0)
            offset=mean(mdsdata(['.BOARD_0' int2str(board) ':CHANNEL_0' num2str(channel,'%0.2d')...
                '[' int2str(tminoffsetindex) ':' int2str(tmaxoffsetindex) ']']))*10/2^15;
        end
    end
end

function [b]=get_chord(board, channel,tbeginindex,tendindex)
if isempty(tbeginindex)
    b=mdsdata(['.BOARD_0' int2str(board) ':CHANNEL_0' num2str(channel,'%0.2d')])*10/2^15;
else
    b=mdsdata(['.BOARD_0' int2str(board) ':CHANNEL_0' num2str(channel,'%0.2d') '[' int2str(tbeginindex) ':' int2str(tendindex) ']'])*10/2^15;
end

function [offset]=get_offset(board, channel,tminoffsetindex,tmaxoffsetindex)
offset=mean(mdsdata(['.BOARD_0' int2str(board) ':CHANNEL_0' num2str(channel,'%0.2d')...
    '[' int2str(tminoffsetindex) ':' int2str(tmaxoffsetindex) ']']))*10/2^15;