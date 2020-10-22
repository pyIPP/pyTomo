function [filt,abso,error]=XraySignals_Strings(det)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  Analyse string contaning the information about the XRAY 
%  diods. Prepares data for XraySignals.m
% 
%  SYNTAX
%   function [filt,abso,error]=XraySignals_Strings(det);
%     
%  INPUT string
%    EXAMPLE det=['Be50,  Si0.7, O1500, N4500 $ Be200']
%
%from Alexei Zabolotsky routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error=0;

ok=~isspace(det);det=det(ok);            % remove blanc spaces  
ind=findstr(det,'$');
if isempty(ind)
   disp('Strings ERROR: No detector found in det string')   
   error=1;
   return
end
if length(ind)>1
   disp('Strings ERROR: More than one detector found in det string')   
   error=1;
   return
end    
   
filtersstr=det(1:ind-1); % take only symbols before "$" sign
detectorstr=det(ind+1:end);

%*************************** Translate detector string ************************************* 

index=[0,findstr(detectorstr,','),length(detectorstr)+1]; 
nabso=length(findstr(detectorstr,','))+1;       % number of absorbers in detector
disp('      ****************  Absorbers: ********************')
for i=1:nabso
   tmp2=detectorstr(index(i)+1:index(i+1)-1);
   abso(i).thick=str2num(tmp2(~isletter(tmp2)));  % take numbers as thickness 
   abso(i).material=tmp2(isletter(tmp2));             % take letters as absorber
   disp(['      *         ',abso(i).material,'  thickness ',int2str(abso(i).thick),'  micron'])
end   

%*************************** Translate filter string ************************************* 

index=[0,findstr(filtersstr,','),length(filtersstr)+1]; 
nfilters=length(findstr(filtersstr,','))+1;       % number of filters
disp('      ****************  Filters: *********************')
for i=1:nfilters
   tmp2=filtersstr(index(i)+1:index(i+1)-1);
   filt(i).thick=str2num(tmp2(~isletter(tmp2)));  % take numbers as thickness 
   filt(i).material=tmp2(isletter(tmp2));             % take letters as absorber
   disp(['      *         ',filt(i).material,'  thickness ',int2str(filt(i).thick),'  micron'])
end   
disp('      ************************************************')

