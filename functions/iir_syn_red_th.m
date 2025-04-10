function ret=iir_syn_red_th(Am,Su,ii,q,eps)

% Am - Theoretically VAR coeficients matrix
% Su - Theoretically noise covariance matrix
% ii - Target process index
% q - Truncation lag for correlations
% eps - pre defined tolerance

%% 1 -- Find MIN -- In this case find synergy 
% Parameters
eIIRm=[]; % Vector to save values of IIR and cIIR 

M=size(Am,1); % Number of Processes
proc_rem_vec=setdiff(1:M,ii); % Remove ii processes from the vector of possibilities to MIN/MAX IIR
comb_zx=nchoosek(proc_rem_vec,2);

% 1st step - Find triplet that MIN the IIR
for i=1:size(comb_zx,1)
    iir_tmp(i)=iir_lin(Am,Su,q,ii,comb_zx(i,1),comb_zx(i,2));
end
[iir_min,imin]=min(iir_tmp);
imin=imin(1);

% Add the two processes that MIN the IIR jointly with the target processe
% ii
in_triplet=[ii comb_zx(imin,1) comb_zx(imin,2)];
% The first two entries of this vector are the first triplet that MIN IIR

if iir_min>=0 % If the first values found is greater than zero no synergy present
    eIIRm=0;
    cm=[];
    in_triplet=[];
else
    % Add value of IIR
    eIIRm=[eIIRm iir_min];

    % Before go for WHILE loop, renew the proc_rem_vect - Possibilities to
    % conditioning
    cset=setdiff(proc_rem_vec,in_triplet);
    exitcrit=0; % Flag for WHILE loop
    cm=[]; % Vector to save indexes of processes

    while exitcrit==0
        IIRtmp=nan*ones(length(cset),1);
        for i=1:length(cset)
            IIRtmp(i)=cIIR_lin(Am,Su,q,in_triplet,[cm cset(i)]);
        end
        [IIRmin,imin]=min(IIRtmp);
        imin=imin(1);

        % If this value is greater than the value before break loop;
        % Otherwise go for surrogate anayisis
        % Also we need to consider when the measure is positive so we have
        % redundancy not synergy
        if IIRmin>eIIRm(end) || IIRmin>=0 % Shift from synergy to redundancy
            break
        end

        eIIRm=[eIIRm IIRmin];
        cm=[cm cset(imin)];
        if abs(IIRmin-eIIRm(end-1))>=eps % cIIR has decreased non-negligibly
            cset(imin)=[]; % Remove the process of cset
            if isempty(cset), exitcrit=1; end
        else % no decrease
            cm(end)=[];
            eIIRm(end)=[];
            exitcrit=1;
        end
    end
    % Struct with values - vector with values of IIR and cIIR, Surrogates and
    % vector with conditioning vector and the initial triplet
end
ret.IIRm=eIIRm;
ret.IIRm_in_triplet=in_triplet;
ret.IIRm_cond_vec=cm;

%% 2 -- Find MAX -- In this case redundancy 
% Parameters
eIIRM=[]; % Vector to save values of IIR 

proc_rem_vec=setdiff(1:M,ii); % Remove ii processes from the vector of possibilities to MIN/MAX IIR

% 1st step - Find triplet that MAX the IIR
for i=1:size(comb_zx,1)
    iir_tmp(i)=iir_lin(Am,Su,q,ii,comb_zx(i,1),comb_zx(i,2));
end
[iir_max,imax]=max(iir_tmp);
imax=imax(1);

% Add the two processes that MAX the IIR jointly with the target processe
% ii
in_triplet=[ii comb_zx(imax,1) comb_zx(imax,2)];
% The first two entries of this vector are the first triplet that MIN IIR

if iir_max<=0
    eIIRM=0;
    cM=[];
    in_triplet=[];
else
    % Add value of IIR for MAX
    eIIRM=[eIIRM iir_max];

    cset=setdiff(proc_rem_vec,in_triplet); % Reset the conditioning vector
    exitcrit=0; % Flag for WHILE loop
    cM=[]; % Vector to store the processes chosen for the conditioning vector
 

    while exitcrit==0
        IIRtmp=nan*ones(length(cset),1);
        for i=1:length(cset)
            IIRtmp(i)=cIIR_lin(Am,Su,q,in_triplet,[cM cset(i)]);
        end
        [IIRmax,imax]=max(IIRtmp);
        imax=imax(1);

        % If this value is less than the value before break loop;
        % Otherwise go for surrogate analysis
        % If IIRmax is negative break ---- this indicate the presence (predominantly) of
        % synergy not redundant
        if IIRmax<eIIRM(end) || IIRmax<=0 % Shift from Redundancy to Synergy
            break
        end

        eIIRM=[eIIRM IIRmax];
        cM=[cM cset(imax)];
        if abs(IIRmax-eIIRM(end-1))>=eps % cIIR has increased non-negligibly (surrogates)
            cset(imax)=[]; % Remove the process of cset
            if isempty(cset), exitcrit=1; end
        else % no decrease
            cM(end)=[];
            eIIRM(end)=[];
            exitcrit=1;
        end
    end
end
ret.IIRM=eIIRM;
ret.IIRM_in_triplet=in_triplet;
ret.IIRM_cond_vec=cM;
end