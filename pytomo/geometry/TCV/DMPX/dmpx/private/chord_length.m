% R,Z ->  flux surface description in cylindrical coordinates
% rc,zc -> coordinates of two extreme points of the chord

function [D]=chord_length(R,Z,rc,zc)
global_p;

[XI,YI]=polyxpoly(R,Z,rc,zc);
i=[1:length(XI)];
D=sqrt((XI(rem(i,2)==0)-XI(rem(i,2)==1)).^2+(YI(rem(i,2)==0)-YI(rem(i,2)==1)).^2);
if isempty(D)
 D=0;
end
