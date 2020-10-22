function rho_calc
global_p
K=(chordZ(1,n_rmin:n_rmax)-chordZ(2,n_rmin:n_rmax))...
./(chordR(1,n_rmin:n_rmax)-chordR(2,n_rmin:n_rmax));
L=chordZ(2,n_rmin:n_rmax)-K.*chordR(2,n_rmin:n_rmax);
rtmp=find(max(psi_1)>=0);
mzpsi=max(psi_1');
ztmp=find(mzpsi>=0);
nminR=rtmp(1);
nmaxR=rtmp(length(rtmp))+1;
[m,n]=max(mzpsi);
[m,n]=min(mzpsi(1:n));
nminZ=n;
nmaxZ=ztmp(length(ztmp))+1;
psi_Rtmp=psi_R(nminR:nmaxR);
psi_Ztmp=psi_Z(nminZ:nmaxZ);
psi_1tmp=psi_1(nminZ:nmaxZ,nminR:nmaxR);
DRI=(psi_Rtmp(2)-psi_Rtmp(1))/20;
RI=psi_Rtmp(1):DRI:psi_Rtmp(length(psi_Rtmp));
DZI=(psi_Ztmp(2)-psi_Ztmp(1))/10;
ZI=psi_Ztmp(1):DZI:psi_Ztmp(length(psi_Ztmp));
FLUXI=interp2(psi_Rtmp,psi_Ztmp,psi_1tmp,RI,ZI','cubic');
for j=1:length(ZI)
 RR(j,:)=(ZI(j)-L)./K;
end
RRR=round((RR-RI(1)+DRI)/DRI);
[fm,fn]=size(FLUXI);
[m,n]=find(RRR>fn);
for i=1:length(m);
 RRR(m(i),n(i))=fn;
end 
for j=1:length(ZI)-40
 xx(j,n_rmin:n_rmax)=FLUXI(j,RRR(j,n_rmin:n_rmax));
end
[mm,nn]=max(xx); 
rho=sqrt(1-mm/max(mm));
[ll,ii]=min(rho);
rho(1:ii)=-rho(1:ii);
for i=n_rmin:n_rmax
 R(i)=RR(nn(i),i);
end
return;
