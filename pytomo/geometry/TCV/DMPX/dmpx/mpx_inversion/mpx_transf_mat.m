%mpx_transf_mat.m calculates the transfer matrix T between soft X-ray
%	emissivity g(rho) and the stationary signal s received by the 
%	MPX chords. The pixels are limited by flux surfaces. The 
%	slit of the DMPX is taken into account.
%
%			s = T*g
%
%Calculate the transfert matrix T between soft X-ray emissivity g [W/m^3] 
%and the power P [W] received by the MPX chords.
%		P=T.g 		P: n_chords*1, 
%					T: n_chords*n_pixel,
%					g: n_pixel*1
%A pixel lies between two flux surfaces (N_FS=n_pixel+1)
%
%	[T]=mpx_transf_mat(sl,sh,geom,I_chords)
%
%
%Inputs:
%	sl		low intersection of the chords with
%			the flux surfaces. Dim: Nfluxsurf x Nchords
%	sh		contains the high intersection of the chords with
%			the flux surfaces. Dim: Nfluxsurf x Nchords
%	geom	structure containing the geometry
%			see mpxdata.m
%	I_chords	chords passing inside the plasma
%
%Outputs:
%	T		transfer matrix
%


function [T]=mpx_transf_mat(sl,sh,I_chords)

%calculates the chord length in each flux surface
l_chord = (sh(:,I_chords)-sl(:,I_chords))';  % dim : n_chords*n_FS
l_chord(isnan(l_chord)) = 0;
T=diff(l_chord,1,2);
