clear all; close all; clc;

% Add path
addpath('../../functions');

% Load Struct of Random Model
load('Model_Struct_13_Nodes_den_0.3.mat');

% Reorganize the Am matrix
Am=[];
for i=1:Data.popt
    Am=[Am Data.ModelDel(:,:,i)];
end 
Su=eye(13);

%% Parameters
M=13; %%% Number of Processes
N=10000; %%% Length of Realization
model_order_p=Data.popt; % Model Order
tau=ones(1,M); % Time Delay 
q=22; % Auto correlation truncation lag
nsurr=100; % number surrogates
nsim=20; % number of simulations
alpha=0.05;
eps=10^(-4);
measures={'iir','rsi','oir'};

% Waitbars
hw1 = waitbar(0,'Measures...');
hw3= waitbar(0,'Target...');

for meas=1:length(measures)
    for ii=1:M

        % Estimate Measures
        ret=greedy_red_syn_th(Am,Su,ii,q,eps,measures{meas});

        if meas==1
            Redundancy_th(meas,ii)=ret.IIRM(end);
            Redundancy_cM_th{meas,ii}=ret.IIRM_cond_vec;
            Redundancy_in_triplet_th{meas,ii}=ret.IIRM_in_triplet;
            Synergy_in_triplet_th{meas,ii}=ret.IIRm_in_triplet;
            Synergy_th(meas,ii)=ret.IIRm(end);
            Synergy_cm_th{meas,ii}=ret.IIRm_cond_vec;
        elseif meas==2
            Redundancy_th(meas,ii)=ret.RSIM(end);
            Redundancy_cM_th{meas,ii}=ret.RSIM_cond_vec;
            Redundancy_in_triplet_th{meas,ii}=ret.RSIM_in_triplet;
            Synergy_in_triplet_th{meas,ii}=ret.RSIm_in_triplet;
            Synergy_th(meas,ii)=ret.RSIm(end);
            Synergy_cm_th{meas,ii}=ret.RSIm_cond_vec;
        else
            Redundancy_th(meas,ii)=ret.OIRM(end);
            Redundancy_cM_th{meas,ii}=ret.OIRM_cond_vec;
            Redundancy_in_triplet_th{meas,ii}=ret.OIRM_in_triplet;
            Synergy_in_triplet_th{meas,ii}=ret.OIRm_in_triplet;
            Synergy_th(meas,ii)=ret.OIRm(end);
            Synergy_cm_th{meas,ii}=ret.OIRm_cond_vec;
        end
        waitbar(ii/M,hw3);

    end
    waitbar(meas/numel(measures),hw1);
end

% Delete Waitbairs
hw1.delete();
hw3.delete();

save Theoretical_Values_13_Nodes_den_0.3.mat

