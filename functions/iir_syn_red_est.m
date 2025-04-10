function ret=iir_syn_red_est(Y,p,ii,q,numsurr,alpha)

% Y - Data (each series in row)
% p - VAR model p
% ii - Target process index
% q - Truncation lag for correlations
% numsurr - number of surrogates to employ
% alpha - Signficance level for non parametric test - Surrogates

%% 1 -- Find MIN -- IN this case find synergy

% Parameters
eIIRm=[]; % Vector to save values of IIR and co

% Estimate Full Model
[eAm,eSu]=lrp_idMVAR(Y,p);

M=size(eAm,1); % Number of Processes
proc_rem_vec=setdiff(1:M,ii); % Remove ii processes from the vector of possibilities to MIN/MAX IIR
comb_zx=nchoosek(proc_rem_vec,2);

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
    eIIRm=0;
    IIRm_Surr=[];
    cm=[];
    in_triplet=[];
elseif iir_min<=IIR_initial_tmps_th(2) && iir_min>=IIR_initial_tmps_th(1) 
    eIIRm=iir_min;
    cm=[];
    IIRm_Surr=[];
else
    % Add value of IIR
    eIIRm=[eIIRm iir_min];


    % Before go for WHILE loop, renew the proc_rem_vect - Possibilities to
    % conditioning
    cset=setdiff(proc_rem_vec,in_triplet);
    exitcrit=0; % Flag for WHILE loop
    cm=[]; % Vector to save indexes of processes
    IIRm_Surr=[];

    while exitcrit==0
        IIRtmp=nan*ones(length(cset),1);
        for i=1:length(cset)
            IIRtmp(i)=cIIR_lin(eAm,eSu,q,in_triplet,[cm cset(i)]);
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

        IIRtmps=nan*ones(numsurr,1);
        for is=1:numsurr
            y=Y(cset(imin),:)';
            ys=surr_iaafft(y);
            Ys=Y; Ys(cset(imin),:)=ys';
            % [p,~,aic,~,~] = lrp_mos_idMVAR(Ys,15);
            [seAm,seSu]=lrp_idMVAR(Ys,p);
            IIRtmps(is)=cIIR_lin(seAm,seSu,q,in_triplet,[cm cset(imin)]);
        end
        IIRtmps_th=prctile(IIRtmps,100*alpha);
        IIRm_Surr=[IIRm_Surr IIRtmps];

        eIIRm=[eIIRm IIRmin];
        cm=[cm cset(imin)];
        if IIRtmps_th>IIRmin % cIIR has decreased non-negligibly (surrogates)
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
ret.IIRm_Surr=IIRm_Surr;
ret.IIRm_in_triplet=in_triplet;
ret.IIRm_cond_vec=cm;

%% 2 -- Find MAX -- In this case 
% Parameters
eIIRM=[]; % Vector to save values of IIR and co

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
IIR_initial_tmps_th=prctile(IIR_initial_tmps,[100*alpha 100*(1-alpha)]);

if iir_max<=0
    eIIRM=0;
    cM=[];
    in_triplet=[];
    IIRM_Surr=[];
elseif iir_max<=IIR_initial_tmps_th(2) && iir_max>=IIR_initial_tmps_th(1) 
    eIIRM=iir_max;
    cM=[];
    IIRM_Surr=[];
else
    % Add value of IIR for MAX
    eIIRM=[eIIRM iir_max];

    cset=setdiff(proc_rem_vec,in_triplet); % Reset the conditioning vector
    exitcrit=0; % Flag for WHILE loop
    cM=[]; % Vector to
    IIRM_Surr=[];
   

    while exitcrit==0
        IIRtmp=nan*ones(length(cset),1);
        for i=1:length(cset)
            IIRtmp(i)=cIIR_lin(eAm,eSu,q,in_triplet,[cM cset(i)]);
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

        IIRtmps=nan*ones(numsurr,1);
        for is=1:numsurr
            y=Y(cset(imax),:)';
            ys=surr_iaafft(y);
            Ys=Y; Ys(cset(imax),:)=ys';
            [seAm,seSu]=lrp_idMVAR(Ys,p);
            IIRtmps(is)=cIIR_lin(seAm,seSu,q,in_triplet,[cM cset(imax)]);
        end
        IIRtmps_th=prctile(IIRtmps,100*(1-alpha));
        IIRM_Surr=[IIRM_Surr IIRtmps];

        eIIRM=[eIIRM IIRmax];
        cM=[cM cset(imax)];
        if IIRtmps_th<IIRmax % cIIR has increased non-negligibly (surrogates)
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
ret.IIRM_Surr=IIRM_Surr;
ret.IIRM_in_triplet=in_triplet;
ret.IIRM_cond_vec=cM;
end