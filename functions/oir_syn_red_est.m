function ret=oir_syn_red_est(Y,p,ii,q,numsurr,alpha)

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
eOIRm=[]; % Vector to save values of IIR and co

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
    eOIRm=0;
    OIRm_Surr=[];
    cm=[];
    in_triplet=[];
elseif iir_min<=IIR_initial_tmps_th(2) && iir_min>=IIR_initial_tmps_th(1) % For now I put this condtion, but we should implement surrogates!!
    eOIRm=iir_min;
    cm=[];
    OIRm_Surr=[];
else
    % Add value of IIR
    eOIRm=[eOIRm iir_min];


    % Before go for WHILE loop, renew the proc_rem_vect - Possibilities to
    % conditioning
    cset=setdiff(1:M,in_triplet);
    exitcrit=0; % Flag for WHILE loop
    cm=in_triplet; % Vector to save indexes of processes
    OIRm_Surr=[];

    while exitcrit==0
        OIRtmp=nan*ones(length(cset),1);
        for i=1:length(cset)
            OIRtmp(i)=OIR_lin(eAm,eSu,q,[cm cset(i)]);
        end
        [OIRmin,imin]=min(OIRtmp);
        imin=imin(1);

        % If this value is greater than the value before break loop;
        % Otherwise go for surrogate anayisis
        % Also we need to consider when the measure is positive so we have
        % redundancy not synergy
        if OIRmin>eOIRm(end) || OIRmin>=0 % Shift from synergy to redundancy
            break
        end

        OIRtmps=nan*ones(numsurr,1);
        for is=1:numsurr
            y=Y(cset(imin),:)';
            ys=surr_iaafft(y);
            Ys=Y; Ys(cset(imin),:)=ys';
            [seAm,seSu]=lrp_idMVAR(Ys,p);
            OIRtmps(is)=OIR_lin(seAm,seSu,q,[cm cset(imin)]);
        end
        OIRtmps_th=prctile(real(OIRtmps),100*alpha);
        OIRm_Surr=[OIRm_Surr OIRtmps];

        eOIRm=[eOIRm OIRmin];
        cm=[cm cset(imin)];
        if OIRtmps_th>OIRmin % cIIR has decreased non-negligibly (surrogates)
            cset(imin)=[]; % Remove the process of cset
            if isempty(cset), exitcrit=1; end
        else % no decrease
            cm(end)=[];
            eOIRm(end)=[];
            exitcrit=1;
        end
    end
end
ret.OIRm=eOIRm;
ret.OIRm_Surr=OIRm_Surr;
ret.OIRm_in_triplet=in_triplet;
ret.OIRm_Final_Vec=cm;

%% 2 -- Find MAX -- In this case 
% Parameters
eOIRM=[]; % Vector to save values of IIR and co

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
    eOIRM=0;
    cM=[];
    in_triplet=[];
    OIRM_Surr=[];
elseif iir_max<=IIR_initial_tmps_th(2) && iir_max>=IIR_initial_tmps_th(1)
    eOIRM=iir_max;
    cM=[];
    OIRM_Surr=[];
else
    % Add value of IIR for MAX
    eOIRM=[eOIRM iir_max];

    cset=setdiff(proc_rem_vec,in_triplet); % Reset the conditioning vector
    exitcrit=0; % Flag for WHILE loop
    cM=in_triplet; % Vector to
    OIRM_Surr=[];

    while exitcrit==0
        OIRtmp=nan*ones(length(cset),1);
        for i=1:length(cset)
            OIRtmp(i)=OIR_lin(eAm,eSu,q,[cM cset(i)]);
        end
        [OIRmax,imax]=max(OIRtmp);
        imax=imax(1);

        % If this value is less than the value before break loop;
        % Otherwise go for surrogate analysis
        % If IIRmax is negative break ---- this indicate the presence (predominantly) of
        % synergy not redundant
        if OIRmax<eOIRM(end) || OIRmax<=0 % Shift from Redundancy to Synergy
            break
        end

        OIRtmps=nan*ones(numsurr,1);
        for is=1:numsurr
            y=Y(cset(imax),:)';
            ys=surr_iaafft(y);
            Ys=Y; Ys(cset(imax),:)=ys';
            [p,~,aic,~,~] = lrp_mos_idMVAR(Ys,15);
            [seAm,seSu]=lrp_idMVAR(Ys,p);
            OIRtmps(is)=OIR_lin(seAm,seSu,q,[cM cset(imax)]);
        end
        OIRtmps_th=prctile(real(OIRtmps),100*(1-alpha));
        OIRM_Surr=[OIRM_Surr OIRtmps];

        eOIRM=[eOIRM OIRmax];
        cM=[cM cset(imax)];
        if OIRtmps_th<OIRmax % OIR has increased non-negligibly (surrogates)
            cset(imax)=[]; % Remove the process of cset
            if isempty(cset), exitcrit=1; end
        else % no decrease
            cM(end)=[];
            eOIRM(end)=[];
            exitcrit=1;
        end
    end
end
ret.OIRM=eOIRM;
ret.OIRM_Surr=OIRM_Surr;
ret.OIRM_in_triplet=in_triplet;
ret.OIRM_Final_Vec=cM;
end