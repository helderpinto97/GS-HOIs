function ret=rsi_syn_red_est(Y,p,ii,q,numsurr,alpha)

% Y - Data (each series in row)
% p - VAR model p
% ii - Target process index
% q - Truncation lag for correlations
% numsurr - number of surrogates to employ
% alpha - Signficance level for non parametric test - Surrogates

% Estimate Full Model
[eAm,eSu]=lrp_idMVAR(Y,p);

M=size(eAm,1); % Number of Processes
proc_rem_vec=setdiff(1:M,ii); % Remove ii processes from the vector of possibilities to MIN/MAX IIR
comb_zx=nchoosek(proc_rem_vec,2);

%% 1 -- Find MIN -- IN this case find synergy so put minus sign

% Parameters
eRSIm=[]; % Vector to save values of IIR and co

% 1st step - Find triplet that MIN the IIR
for i=1:size(comb_zx,1)
    iir_tmp(i)=iir_lin(eAm,eSu,q,ii,comb_zx(i,1),comb_zx(i,2));
end
[iir_min,imin]=min(iir_tmp);
imin=imin(1);

% Add the two processes that MIN the IIR jointly with the target processe
% ii
in_triplet=[ii comb_zx(imin,1) comb_zx(imin,2)];
% The first two entries of this vector are the first triplet that MIN IIR

% Create surrogates to test significance of iir initial value
IIR_initial_tmps=nan*ones(numsurr,1);
for is=1:numsurr
    y=Y(comb_zx(imin,1),:)';
    ys=surr_iaafft(y);
    Ys=Y; Ys(comb_zx(imin,1),:)=ys';
    [se_in_Am,se_in_Su]=lrp_idMVAR(Ys,p);
    IIR_initial_tmps(is)=iir_lin(se_in_Am,se_in_Su,q,ii,comb_zx(imin,1),comb_zx(imin,2));
end
IIR_initial_tmps_th=prctile(real(IIR_initial_tmps),[100*alpha 100*(1-alpha)]);

if iir_min>=0
    eRSIm=0;
    RSIm_Surr=[];
    cm=[];
    in_triplet=[];
elseif iir_min<=IIR_initial_tmps_th(2) && iir_min>=IIR_initial_tmps_th(1) % For now I put this condtion, but we should implement surrogates!!
    eRSIm=iir_min;
    cm=[];
    RSIm_Surr=[];
else
    % Add value of IIR
    eRSIm=[eRSIm iir_min];


    % Before go for WHILE loop, renew the proc_rem_vect - Possibilities to
    % conditioning
    cset=setdiff(1:M,in_triplet);
    exitcrit=0; % Flag for WHILE loop
    cm=in_triplet; % Vector to save indexes of processes
    RSIm_Surr=[];

    while exitcrit==0
        RSItmp=nan*ones(length(cset),1);
        for i=1:length(cset)
            RSItmp(i)=RSI_lin(eAm,eSu,q,[cm cset(i)]);
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

        RSItmps=nan*ones(numsurr,1);
        for is=1:numsurr
            y=Y(cset(imin),:)';
            ys=surr_iaafft(y);
            Ys=Y; Ys(cset(imin),:)=ys';
            [seAm,seSu]=lrp_idMVAR(Ys,p);
            RSItmps(is)=RSI_lin(seAm,seSu,q,[cm cset(imin)]);
        end
        RSItmps_th=prctile(RSItmps,100*alpha);
        RSIm_Surr=[RSIm_Surr RSItmps];

        eRSIm=[eRSIm RSImin];
        cm=[cm cset(imin)];
        if RSItmps_th>RSImin % cIIR has decreased non-negligibly (surrogates)
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
ret.RSIm_Surr=RSIm_Surr;
ret.RSIm_in_triplet=in_triplet;
ret.RSIm_Final_Vec=cm;

%% 2 -- Find MAX -- In this case 
% Parameters
eRSIM=[]; % Vector to save values of IIR and co

proc_rem_vec=setdiff(1:M,ii); % Remove ii processes from the vector of possibilities to MIN/MAX IIR

% 1st step - Find triplet that MAX the IIR
for i=1:size(comb_zx,1)
    iir_tmp(i)=iir_lin(eAm,eSu,q,ii,comb_zx(i,1),comb_zx(i,2));
end
[iir_max,imax]=max(iir_tmp);
imax=imax(1);

% Add the two processes that MAX the IIR jointly with the target processe
% ii
in_triplet=[ii comb_zx(imax,1) comb_zx(imax,2)];
% The first two entries of this vector are the first triplet that MIN IIR

% Do surrogates to test the significance of iir_max value. If not
% statistically significant dont proceed with the surrogates
IIR_initial_tmps=nan*ones(numsurr,1);
for is=1:numsurr
    y=Y(comb_zx(imax,1),:)';
    ys=surr_iaafft(y);
    Ys=Y; Ys(comb_zx(imax,1),:)=ys';
    [se_in_Am,se_in_Su]=lrp_idMVAR(Ys,p);
    IIR_initial_tmps(is)=iir_lin(se_in_Am,se_in_Su,q,ii,comb_zx(imax,1),comb_zx(imax,2));
end
IIR_initial_tmps_th=prctile(real(IIR_initial_tmps),[100*alpha 100*(1-alpha)]);


if iir_max<=0
    eRSIM=0;
    cM=[];
    in_triplet=[];
    RSIM_Surr=[];
elseif iir_max<=IIR_initial_tmps_th(2) && iir_max>=IIR_initial_tmps_th(1)
    eRSIM=iir_max;
    cM=[];
    RSIM_Surr=[];
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
            RSItmp(i)=RSI_lin(eAm,eSu,q,[cM cset(i)]);
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

        RSItmps=nan*ones(numsurr,1);
        for is=1:numsurr
            y=Y(cset(imax),:)';
            ys=surr_iaafft(y);
            Ys=Y; Ys(cset(imax),:)=ys';
            [seAm,seSu]=lrp_idMVAR(Ys,p);
            RSItmps(is)=RSI_lin(seAm,seSu,q,[cM cset(imax)]);
        end
        RSItmps_th=prctile(real(RSItmps),100*(1-alpha));
        RSIM_Surr=[RSIM_Surr RSItmps];

        eRSIM=[eRSIM RSImax];
        cM=[cM cset(imax)];
        if RSItmps_th<RSImax % cIIR has increased non-negligibly (surrogates)
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
ret.RSIM_Surr=RSIM_Surr;
ret.RSIM_in_triplet=in_triplet;
ret.RSIM_Final_Vec=cM;
end