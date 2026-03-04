function dfoutput = dmpxfir_snapshot( shots, pps )
%DMPXFIR_SNAPSHOT dfoutput = dmpxfir_snapshot( shots, pps )
%   For all given shots this function prints a DMPX/FIR time trace. The
%   resolution of the saved PNG file is defined by the variable pps
%   (Pixels per second). The default pps is 880 (FullHD width for ~2.2s).
%
%   J. Kamleitner, CRPP, EPFL, Apr 2012

if(nargin<2)
    pps=880;
end

fprintf('\nDMPX+FIR snapshot program started.\n');


for i=1:length(shots)
    %----------------------------------------------------------------------
    % get data
    %----------------------------------------------------------------------
    fprintf('- shot %d:\n  . reading data ...\n',shots(i));
    tread=tic;
    try
        dmpxdata=mpxdata(shots(i),'s','freq',1,'chords',32);
        td=dmpxdata.signal.dim{1};
    catch ex
        fprintf(ex.message);
        td=[NaN NaN];
    end
    mdsopen(shots(i))
    try
        xtedata=tdi('\results::te_x_a');
        tx=xtedata.dim{1};
    catch ex
        fprintf(ex.message);
        tx=[NaN NaN];
    end
    try
        firdata=tdi('\results::fir:n_average');
        tf=firdata.dim{1};
    catch ex
        fprintf(ex.message);
        tf=[NaN NaN];
    end
    tdiag=[min([min(td) min(tx) min(tf)]) max([max(td) max(tx) max(tf)])];
    if(isnan(tdiag(1)))
        fprintf('  . no data for this shot, continuing with next one ...\n');
        continue;
    end
    try % get liuqe times
        temp=tdi('TCV_EQ("psi")');    % times liuqe all (all available liuqe points in time)
        tla=temp.dim{3};
    catch ex
        fprintf(ex.message);
        tla=[NaN NaN];
    end
    mdsclose;
    tl=[min(tla)-0.02 max(tla)+0.05];
    t=[max(tdiag(1),tl(1)), min(tdiag(2),tl(2))];
    if(diff(t)<0.02)
        fprintf('  . Time interval [%g,%g]s too short, cancelling this shot!\n',t(1),t(2));
        continue;
    end
    % select data in time window
    if(isfinite(td(1)))
        dmpxdat=dmpxdata.signal.data( (td>t(1)) & (td<t(2)) , : );
        td=td( (td>t(1)) & (td<t(2)) , : );
        if(numel(td)<2)
            td=[NaN NaN];
        end
        dmpxmax=max(dmpxdat);
    end
    if(isfinite(tx(1)))
        xtedat=nanmean(xtedata.data( (tx>t(1)) & (tx<t(2)) , : ),2);
        xtestd=nanstd(xtedata.data( (tx>t(1)) & (tx<t(2)) , : ),0,2);
        tx=tx( (tx>t(1)) & (tx<t(2)) , : );
        if(numel(tx)<2 || max(xtestd)==0)
            tx=[NaN NaN];
        end
        xtemax=max(xtedat);
    end
    if(isfinite(tf(1)))
        firdat=firdata.data( (tf>t(1)) & (tf<t(2)) , : );
        tf=tf( (tf>t(1)) & (tf<t(2)) , : );
        if(numel(tf)<2)
            tf=[NaN NaN];
        end
        firmaxlimit=2e20;
        firmax=min(max(firdat),firmaxlimit);
    end
    fprintf('  . data read within %gs.\n',toc(tread));
    
    for j=1:length(pps)
        %----------------------------------------------------------------------
        % init figure
        %----------------------------------------------------------------------
        fprintf('  . initalizing figure.\n');
        [ifig,sf]=snap_simple_initfig('timetrace',t,pps(j),150);

        %----------------------------------------------------------------------
        % plot and save figure
        %----------------------------------------------------------------------
        fprintf('  . plotting figure.\n');
        
        hold on;
        if(isfinite(tx(1)))
            plot(tx, xtedat/xtemax,'m-','LineWidth',2/sf);
            plot(tx, xtestd/xtemax,'m:','LineWidth',1/sf);
        end
        if(isfinite(td(1)))
            plot(td, dmpxdat/dmpxmax,'k-','LineWidth',2/sf);
        end
        if(isfinite(tf(1)))
            plot(tf, firdat/firmax,'b-','LineWidth',2/sf);
        end
        if(isfinite(td(1)))
            text(t(1)+0.02*diff(t),0.88,sprintf('DMPX sig_{ch32}<=%.2g',dmpxmax),'Color','k','FontSize',16/sf,'background','w');
        end
        if(isfinite(tf(1)))
            text(t(1)+0.02*diff(t),0.7,sprintf('FIR n_{e,av,19}<=%3.1f',firmax*1e-19),'Color','b','FontSize',16/sf,'background','w');
        end
        if(isfinite(tx(1)))
            text(t(1)+0.02*diff(t),0.52,sprintf('X-T_{e,eV}<=%4.0f',xtemax),'Color','m','FontSize',16/sf,'background','w');
        end
        hold off;
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        xlim(t);
        ylim([0 1]);

        snap_simple_savefig(ifig,'base/ttFDX/ttFIRDMPX',shots(i),t,pps(j),'dmpxfir_snapshot');

    end
    
end

dfoutput=0;

end

