function init12w
global_p;

fig1=figure('units','normalized','position',[.45 .4 .5 .5],... 
	'resize','on','tag','upperleft',... 
	'visible','off');

figure(fig1);
f1_1=subplot('Position',[.06 .55 .6 .4]);
set(gca,'XLim',[min(t) max(t)]);
f1_2=subplot('Position',[.06 .07 .6 .4]);
f1_3=subplot('Position',[.7 .08 .3 .87]);

 h_x_1 = uicontrol(...
  'BackgroundColor',[ .8 .8 .8 ],...
		'CallBack','',... 
  'Units','normalized','Position',[ .1 .47 .2 .04 ],... 
		'HorizontalAlignment','left',... 
		'String','X1',... 
  'Style','text');
   % pointer Y
 h_y_1 = uicontrol(...
  'BackgroundColor',[ .8 .8 .8 ],...
		'CallBack','',... 
  'Units','normalized','Position',[ .1 .001 .2 .04 ],... 
		'HorizontalAlignment','left',... 
		'String','X2',... 
  'Style','text');
  

		 
	x_auto_1 = uimenu('Label',' auto_X',...                    
		'CallBack',['set(gca,''XLimMode'',''auto'');',...
       'set(f1_1,''XLimMode'',''auto'');' ]); 
                         
	y_auto_1 = uimenu('Label',' auto_Y',...                    
		'CallBack',['set(gca,''YLimMode'',''auto'');',...
       'set(f1_1,''YLimMode'',''auto'');']); 
                         
	zoom_1 = uimenu('Label',' zoom','ForegroundColor',[zoo_m zoo_m zoo_m],... 	                   
		'CallBack',['zoo_m=abs(zoo_m-1);','zoom;',...
    'set(zoom_1,''ForegroundColor'',[zoo_m zoo_m zoo_m]);']);                     

	hold_1 = uimenu('Label',' hold_1',...
  'ForegroundColor',[hold_w(1) hold_w(1) hold_w(1)],...                    
  'CallBack',['hold_w(1)=abs(hold_w(1)-1);',...
	 'set(hold_1,''ForegroundColor'',[hold_w(1) hold_w(1) hold_w(1)]);',...
 	'if hold_w(1)~=1;',...
     'n_r=n_r(1);',...
     'plotter2(0);',...
  'end;']);
 hold_2 = uimenu('Label',' hold_2',...
  'ForegroundColor',[hold_w(2) hold_w(2) hold_w(2)],...                    
  'CallBack',['hold_w(2)=abs(hold_w(2)-1);',...
  'set(hold_2,''ForegroundColor'',[hold_w(2) hold_w(2) hold_w(2)]);',...
		'if hold_w(2)~=1;',...
     'n_t=n_t(1);',... 
     'plotter2(0);',...
  'end;']);
 xval_2 = uimenu('Label',' X_2unit');                      
		 rho_2 = uimenu(xval_2,'Label',' rho',...                    
			 'CallBack',['if hold_w(1)==0;',...
                 'if hold_w(2)==0;',...
                   'xlab_2=1;','n_t_last=0;','plotter2(0);',...
                 'end;',...  
                'end']);
		 R_2 = uimenu(xval_2,'Label',' R',...                    
			 'CallBack',['if hold_w(1)==0;',...
                 'if hold_w(2)==0;',...
                  'xlab_2=2;','n_t_last=0;','plotter2(0);',...
                 'end;',...
                'end']);
		 chan_2 = uimenu(xval_2,'Label',' N_chan',...                    
			 'CallBack',['if hold_w(1)==0;',...
                 'if hold_w(2)==0;',...
                  'xlab_2=3;','plotter2(0);',...
                 'end;',...
                'end']);               	 
 swsp_3 = uimenu('Label',' sw_3 ',...                    
  'CallBack',['sws_3=abs(sws_3-1);',...
		  'if sws_3~=1;',...
      'f1_1=subplot(''Position'',[.06 .55 .6 .4]);',...
      'f1_2=subplot(''Position'',[.06 .07 .6 .4]);',...
      'f1_3=subplot(''Position'',[.7 .08 .3 .87]);',...
      'visio;',...
      'plotter2(0);',...
    'else;',...
      'f1_1=subplot(''Position'',[.06 .55 .9 .4]);',...
      'f1_2=subplot(''Position'',[.06 .07 .9 .4]);',...
      'plotter2(0);',... 
    'end;']);   
reread = uimenu('Label','reread');
	rer_win = uimenu(reread,'Label','window',...                    
		'CallBack',['global_p; subplot(f1_1);tmp=get(gca,''XLim'');',...
              't1=tmp(1); t2=tmp(2);',...
              'y=[]; t=[]; pack;',...
              'load_mpx_data;',...              
              'plotter2(0)']);
% rer_shot = uimenu(reread,'Label','shot',...
%  'CallBack',['shot=input(''Enter shot number: '');',...
%              't1=[]; psi=[]; y=[]; t=[]; pack;',...   
%              'if shot<18093;',...
%                't1=input(''Enter t1: '');',...
%                't2=input(''Enter t2: '');',...
%              'end;',...
%              'load_mpx_data;',...
%              'plotter2(0)']);                                          
	quit_1 = uimenu('Label','Exit',...                    
		'CallBack','clear all; close(get(0,''Children''));');  

set(gcf,'visible','on'); 

set(gcf,'WindowButtonDownFcn',...
   ['if zoo_m==0;',...
     'if round(get(gcf,''CurrentAxes''))==round(f1_1);',...
      'tmp=get(gca,''CurrentPoint'');',...
      '[tmp,n_t1]=min(abs(t-tmp(1,1)));',...
      'if hold_w(2)==0;',...
        'n_t=n_t1;',...
      'else;',...
        'n_t=[n_t n_t1];',...
       'end;',... 
     'else;',...
      'tmp=get(gca,''CurrentPoint'');',...
      '[tmp,n_r1]=min(abs(r-tmp(1,1)));',...
      'if hold_w(1)==0;',...
        'n_r=n_r1;',...
      'else;',...
        'n_r=[n_r n_r1];',...
       'end;',... 
     'end;',... 
     'plotter2(0);',...
     'end;']);


axis('auto');
drawnow
