function ret=rsi_syn_red_th(Am,Su,ii,q,eps)

% Am - Theoretical VAR coeficients matrix 
% Su - Theoretical noise covariance matrix
% ii - Target process index
% q - Truncation lag for correlations
% eps - pre defined tolerance

M=size(Am,1); % Number of Processes
proc_rem_vec=setdiff(1:M,ii); % Remove ii processes from the vector of possibilities to MIN/MAX IIR
comb_zx=nchoosek(proc_rem_vec,2);

%% 1 -- Find MIN -- IN this case find synergy

% Parameters
eRSIm=[]; % Vector to save values of IIR and co

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

if iir_min>=0
    eRSIm=0;
    RSIm_Surr=[];
    cm=[];
    in_triplet=[];
else
    % Add value of IIR
    eRSIm=[eRSIm iir_min];


    % Before go for WHILE loop, renew the proc_rem_vect - Possibilities to
    % conditioning
    cset=setdiff(1:M,in_triplet);
    exitcrit=0; % Flag for WHILE loop
    cm=in_triplet; % Vector to save indexes of processes

    while exitcrit==0
        RSItmp=nan*ones(length(cset),1);
        for i=1:length(cset)
            RSItmp(i)=RSI_lin(Am,Su,q,[cm cset(i)]);
        end
        [RSImin,imin]=min(RSItmp);
        imin=imin(1);

        % If this value is greater than the value before break loop;
        % Otherwise go for surrogate anayisis
        % Also we need to consider when the measure is positive so we have
        % redundancy not synergy
        if RSImin>eRSIm(end) || RSImin>=0 % Shift from synergy to redundancy
            break
        end

        eRSIm=[eRSIm RSImin];
        cm=[cm cset(imin)];
        if abs(RSImin-eRSIm(end-1))>=eps % RSI has decreased non-negligibly (surrogates)
            cset(imin)=[]; % Remove the process of cset
            if isempty(cset), exitcrit=1; end
        else % no decrease
            cm(end)=[];
            eRSIm(end)=[];
            exitcrit=1;
        end
    end
end
ret.RSIm=eRSIm;
ret.RSIm_in_triplet=in_triplet;
ret.RSIm_cond_vec=cm;

%% 2 -- Find MAX -- In this case redundancy 
% Parameters
eRSIM=[]; % Vector to save values of IIR and co

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

if iir_max<=0
    eRSIM=0;
    cM=[];
    in_triplet=[];
else
    % Add value of IIR for MAX
    eRSIM=[eRSIM iir_max];

    cset=setdiff(proc_rem_vec,in_triplet); % Reset the conditioning vector
    exitcrit=0; % Flag for WHILE loop
    cM=in_triplet; % Vector to
    RSIM_Surr=[];

    while exitcrit==0
        RSItmp=nan*ones(length(cset),1);
        for i=1:length(cset)
            RSItmp(i)=RSI_lin(Am,Su,q,[cM cset(i)]);
        end
        [RSImax,imax]=max(RSItmp);
        imax=imax(1);

        % If this value is less than the value before break loop;
        % Otherwise go for surrogate analysis
        % If IIRmax is negative break ---- this indicate the presence (predominantly) of
        % synergy not redundant
        if RSImax<eRSIM(end) || RSImax<=0 % Shift from Redundancy to Synergy
            break
        end

        eRSIM=[eRSIM RSImax];
        cM=[cM cset(imax)];
        if abs(RSImax-eRSIM(end-1))>=eps % cIIR has increased non-negligibly (surrogates)
            cset(imax)=[]; % Remove the process of cset
            if isempty(cset), exitcrit=1; end
        else % no decrease
            cM(end)=[];
            eRSIM(end)=[];
            exitcrit=1;
        end
    end
end
ret.RSIM=eRSIM;
ret.RSIM_in_triplet=in_triplet;
ret.RSIM_cond_vec=cM;
end