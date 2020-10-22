%addpath /home/kamleitn/matlab/MatlabMPI/src
%addpath /home/duval/matlab7_5_2010c

% Initialize MPI.
MPI_Init;

% Create communicator.
comm = MPI_COMM_WORLD;

% Get size and rank.
comm_size = MPI_Comm_size(comm);
my_rank = MPI_Comm_rank(comm);


if(my_rank==0)
    % broadcast gti parameter sets
    fprintf('- Broadcasting parameters.\n');
    tbc=tic;
    MPI_Bcast( 0, 1, comm, shot, ichan, DTNE1, DTNE2, ch, echant);
    fprintf('- Parameters broadcasted within %g s.\n',toc(tbc));
    %if(shot<43083)
    %    mdsdisconnect;
    %    fprintf('- connecting to TCV1.\n');
    %    mdsconnect('tcv1');
    %    mdsopen(sprintf('tcv_shot:%d',shot));
    %end
else
    % receive gti parameter sets
    fprintf('- Waiting for parameters.\n');
    trcv=tic;
    [shot, ichan, DTNE1, DTNE2, ch, echant]=MPI_Recv(0,1,comm);
    fprintf('- Parameters received within %g s.\n',toc(trcv));
    if(shot<43083)
        fprintf('- connecting to TCV1.\n');
        mdsconnect('tcv1');
    end
    mdsopen(shot);
end

k=0;
for i=1:length(ichan)
    %fprintf('test\n');
    if(my_rank==mod(i-1,comm_size))
        k=k+1;
        fprintf('- Reading ch # %d.\n',i);
        if(ichan(i)<=32)
            ttmpdata(:,k)=mdsdata([DTNE1 ch(ichan(i),:) echant]);
        else
            ttmpdata(:,k)=mdsdata([DTNE2 ch(ichan(i)-32,:) echant]);
        end
    end
end

if(my_rank>0)
    mdsclose();
    mdsdisconnect;
    % send results to root
    MPI_Send(0,1000+my_rank,comm,ttmpdata)
    fprintf('- Sent data to root.\n');
else
    for ir=0:(comm_size-1)
        k=0;
        if(ir>0)
            % receive computation results from non-root processes
            ttmpdata=MPI_Recv(ir,1000+ir,comm);
            fprintf('- Received data from process # %d.\n',ir);
        end
        for i=1:length(ichan)        
            if(ir==mod(i-1,comm_size))
                k=k+1;
                tmpdata(:,ichan(i))=ttmpdata(:,k);
                fprintf('- Stored ch # %d.\n',i);
            end
        end
    end
end


% Finalize Matlab MPI.
MPI_Finalize;
if (my_rank ~= MatMPI_Host_rank(comm))
  exit;
end