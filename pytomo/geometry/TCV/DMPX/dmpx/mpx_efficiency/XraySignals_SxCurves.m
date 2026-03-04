function [tr,error]=XraySignals_SxCurves(abso,thick,e)

error=0;tr=[];
thick=thick/1e4;  % from microns to cm

abso=upper(abso);
filelocate=which('XraySignals_SxCurves.m');
ind=findstr(filelocate,'XraySignals_SxCurves.m');
if or(strcmp(computer,'IBM_RS'),strcmp(computer,'GLNXA64'))
     directory=[filelocate(1:ind-1) 'DATA/'];
else
      directory=[filelocate(1:ind-1) 'DATA\'];
end     
try
		if 	findstr(abso,'BE')
            out=load([directory 'be.dat']);	
			ev=(out(:,1));
			murho=(out(:,2));     % mass attenuation coefficient (cm2 g-1)
            rho=1.845;            % Nominal density:    (g cm-3)
            tr=exp(-(murho.*rho.*thick));
            
		elseif	findstr(abso,'AL')
         	out=load([directory 'al.dat']);	
			ev=(out(:,1));
			murho=(out(:,2));
            rho=2.6941;
			tr=exp(-(murho.*rho.*thick));
            
		elseif	findstr(abso,'SI')
         	out=load([directory 'si.dat']);	
			ev=(out(:,1));
			murho=(out(:,2));
            rho=2.3200;
            tr=exp(-(murho.*rho.*thick));            
            
		elseif 	findstr(abso,'HE')
         	out=load([directory 'he.dat']);		
			ev=(out(:,1));
            rho=1.6640E-04;
            murho=(out(:,2));
            tr=exp(-(murho.*rho.*thick)); 
            
		elseif 	strcmp(abso,'KRCH')
			out=load([directory 'kr.dat']);		
			ev=(out(:,1));
            rhokr=3.4840E-03;
            murhokr=(out(:,2));
            murhokr=interp1(ev*1000,murhokr,e);

            out=load([directory 'ch4.dat']);
			ev=(out(:,1));
            rhoch4=0.424;
            murhoch4=(out(:,2));
            murhoch4=interp1(ev*1000,murhoch4,e);
 
            %Mixture 90% Kr and 10% CH4 
			tr=exp(-(murhoch4.*rhoch4.*thick*0.1+murhokr.*rhokr.*thick*0.9));     %only KR generates photoelectron detected by MPX wires
 	
		elseif 	strcmp(abso,'KRCHABS')
			out=load([directory 'kr.dat']);		
			ev=(out(:,1));
            rhokr=3.4840E-03;
            murhokr=(out(:,2));
            murhokr=interp1(ev*1000,murhokr,e);

            %Mixture 90% Kr and 10% CH4 
			tr=exp(-(murhokr.*rhokr.*thick*0.9));     
            
        elseif 	strcmp(abso,'ARCH')
         	out=load([directory 'ar.dat']);		
			ev=(out(:,1));
            rhoar=3.4840E-03;
            murhoar=(out(:,2));
            murhoar=interp1(ev*1000,murhoar,e);

            out=load([directory 'ch4.dat']);
			ev=(out(:,1));
            rhoch4=0.424;
            murhoch4=(out(:,2));
            murhoch4=interp1(ev*1000,murhoch4,e);
 
            %Mixture 90% Ar and 10% CH4 
			tr=exp(-(murhoch4.*rhoch4.*thick*0.1+murhoar.*rhoar.*thick*0.9));                
 
        elseif 	strcmp(abso,'ARCHABS')
         	out=load([directory 'ar.dat']);		
			ev=(out(:,1));
            rhoar=3.4840E-03;
            murhoar=(out(:,2));
            murhoar=interp1(ev*1000,murhoar,e);

            %Mixture 90% Ar and 10% CH4 
			tr=exp(-(murhoar.*rhoar.*thick*0.9));                
 		
		elseif 	strcmp(abso,'XECH')
         	out=load([directory 'xe.dat']);		
			ev=(out(:,1));
            rhoxe=5.4580E-03;
            murhoxe=(out(:,2));
            murhoxe=interp1(ev*1000,murhoxe,e);

            out=load([directory 'ch4.dat']);
			ev=(out(:,1));
            rhoch4=0.424;
            murhoch4=(out(:,2));
            murhoch4=interp1(ev*1000,murhoch4,e);
 
            %Mixture 90% Xe and 10% CH4 
			tr=exp(-(murhoch4.*rhoch4.*thick*0.1+murhoxe.*rhoxe.*thick*0.9));     
            
		elseif 	strcmp(abso,'XECHABS')
         	out=load([directory 'xe.dat']);		
			ev=(out(:,1));
            rhoxe=5.4580E-03;
            murhoxe=(out(:,2));
            murhoxe=interp1(ev*1000,murhoxe,e);

            %Mixture 90% Xe and 10% CH4 
			tr=exp(-(murhoxe.*rhoxe.*thick*0.9));     
  
  		elseif 	findstr(abso,'AIR')          
            out=load([directory 'n.dat']);		
			ev=(out(:,1));
            rhon=1.1650E-03;
            murhon=(out(:,2));
            murhon=interp1(ev*1000,murhon,e);

            out=load([directory 'o.dat']);		
			ev=(out(:,1));
            rhoo=1.3310E-03;
            murhoo=(out(:,2));
            murhoo=interp1(ev*1000,murhoo,e);
 
            out=load([directory 'c.dat']);		
			ev=(out(:,1));
            rhoc=2.2600;
            murhoc=(out(:,2));
            murhoc=interp1(ev*1000,murhoc,e);

         	out=load([directory 'ar.dat']);		
			ev=(out(:,1));
            rhoar=3.4840E-03;
            murhoar=(out(:,2));
            murhoar=interp1(ev*1000,murhoar,e);
            
            %Mixture 78.1% N, 21% O, 0.015% C, 0.47% Ar 
			tr=exp(-(murhon.*rhon.*thick*0.781+murhoo.*rhoo.*thick*0.21+...
                     murhoc.*rhoc.*thick*0.00015+murhoar.*rhoar.*thick*0.0047));     
            
        else
            disp(['SxCurves ERROR: Data for ' abso ' does not exist'])
            error=1;
            return
		end 
catch
        disp(['SxCurves ERROR: cannot load data file for ' abso]);
        error=1;
end    

if length(tr)~=length(e)
     tr=interp1(ev*1000,tr,e);   
end            
