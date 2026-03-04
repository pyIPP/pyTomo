%Test whether the mpx fast data acquisition is available in the MDS tree
%or not (test only the first channel of each card)
%IF NOT, ask for it to the DMPX responsible officer (Loic CURCHOD 10.08.2007)
%
%[test]=isfastmpx(shot)
%
%input: 	shot: 	list of shots
%output: 	test: 	1->fast acquisition is available
%			0-> only slow acquisition (20kHz) is available in the MDS tree
%
% Note: only the first channel of the first card is tested

% Written by Yann Camenen.
% Modified by Loic Curchod, 03.03.2009. To take a vector as input.


function [test]=isfastmpx(shot)

test=NaN;
for ii=1:length(shot),
	if shot(ii) < 20030,
		disp('This program works only for shot>=20030');
	elseif shot(ii) > 34987
		test(ii) = 1;	% Since 34988, all shots are stored with 200 kHz acquisition
	else
		test(ii)=(tt(shot(ii)));	% Il faut passer par une sous fonction pour reinitialiser TDI,
						% sinon au bout d'un moment, il plante....
	end
end

%------------------------------
function test_fast=tt(s)

mdsopen(s)
if (s<24087)|(s>24725)
	l1=mdsdata('node_size("\\ATLAS::DT100_NORTHEAST_001:SELECTED:CHANNEL_001")'); 
	l2=mdsdata('node_size("\\ATLAS::DT100_NORTHEAST_002:SELECTED:CHANNEL_001")');
	test_fast=(l1>0)&(l2>0);
	if s>=26555
		l3=mdsdata('node_size("\\ATLAS::DT100_NORTHEAST_003:SELECTED:CHANNEL_001")');
		test_fast=test_fast&(l3>0);
	end
else
	l1=mdsdata('node_size("\\ATLAS::DT100_SOUTHWEST_001:SELECTED:CHANNEL_001")');
	l2=mdsdata('node_size("\\ATLAS::DT100_SOUTHWEST_002:SELECTED:CHANNEL_001")');
	test_fast=(l1>0)&(l2>0);
end
