% Small routine to help the DJ setting the gains for the DMPX
%
%	DMPX_HV_for_DJ(mode)
%
% mode:	0 -> give the HV value used for a given shot
%       1 -> calculate the HV value to obtain the desired signal level
%      >1 -> use the mode directly as shot number

function []=DMPX_HV_for_DJ(mode)

if(mode>1)
    fprintf('\nUsing mode directly as shot number.\n');
    shot=mode;
    mode=0;
elseif(mode==0)
    shot=input('Shot number: ');
end

switch mode
case 0
    res=mpxdata(shot,'v');
    disp(['Shot: ' num2str(shot)])
    disp(['Top detector: ' num2str(res.top.voltage) 'V'])	
    disp(['Bottom detector: ' num2str(res.bot.voltage) 'V'])	
case 1
    Vref=input('Present high voltage value: ');
    Sref=input('Present signal level in the DJ scope window: ');
    Sout=input('Desired signal level (between 0 and 10V): ');
    Vout=10*round(200*log(Sout*exp(0.005*Vref)/Sref)/10);
    if(Vout>2500)
        Vout=2500;
    end
    disp(['Set the high voltage value to ' num2str(Vout) 'V'])	

end

end
