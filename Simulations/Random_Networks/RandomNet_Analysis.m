clear all; close all; clc;

% Add path
addpath('../../functions');

% Load Struct of Random Model
load('Model_Struct_13_Nodes_den_0.1.mat');

% Reorganize the Am matrix
Am=[];
for i=1:Data.popt
    Am=[Am Data.ModelDel(:,:,i)];
end 

%% Parameters
M=13; %%% Number of Processes
N=10000; %%% Length of Realization
model_order_p=Data.popt; % Model Order
tau=ones(1,M); % Time Delay 
q=22; % Auto correlation truncation lag
nsurr=100; % number surrogates
nsim=20; % number of simulations
alpha=0.05;
measures={'iir','rsi','oir'};

% Waitbars
hw1 = waitbar(0,'Measures...');
hw2= waitbar(0,'Simulation Loop...');
hw3= waitbar(0,'Target...');

for meas=1:length(measures)
    for isim=1:nsim
        U = mvnrnd(zeros(1,M),Data.Sw_nolag0,N);
        Y=var_filter(Am,U);

        % Normalize the Data
        Y=normalize(Y);
        % Series on row
        Y=Y';
        for ii=1:M

            % Estimate Measures
            ret=greedy_red_syn_est(Y,Data.popt,ii,q,nsurr,alpha,measures{meas});

            if meas==1
                Redundancy_est(isim,ii,meas)=ret.IIRM(end);
                Redundancy_cM{isim,ii,meas}=ret.IIRM_cond_vec;
                Redundancy_in_triplet{isim,ii,meas}=ret.IIRM_in_triplet;
                Synergy_in_triplet{isim,ii,meas}=ret.IIRm_in_triplet;
                Synergy_est(isim,ii,meas)=ret.IIRm(end);
                Synergy_cm{isim,ii,meas}=ret.IIRm_cond_vec;
            elseif meas==2
                Redundancy_est(isim,ii,meas)=ret.RSIM(end);
                Redundancy_cM{isim,ii,meas}=ret.RSIM_Final_Vec;
                Redundancy_in_triplet{isim,ii,meas}=ret.RSIM_in_triplet;
                Synergy_in_triplet{isim,ii,meas}=ret.RSIm_in_triplet;
                Synergy_est(isim,ii,meas)=ret.RSIm(end);
                Synergy_cm{isim,ii,meas}=ret.RSIm_Final_Vec;
            else
                Redundancy_est(isim,ii,meas)=ret.OIRM(end);
                Redundancy_cM{isim,ii,meas}=ret.OIRM_Final_Vec;
                Redundancy_in_triplet{isim,ii,meas}=ret.OIRM_in_triplet;
                Synergy_in_triplet{isim,ii,meas}=ret.OIRm_in_triplet;
                Synergy_est(isim,ii,meas)=ret.OIRm(end);
                Synergy_cm{isim,ii,meas}=ret.OIRm_Final_Vec;
            end
            waitbar(ii/M,hw3);
        end
        waitbar(isim/nsim,hw2);
    end
    waitbar(meas/numel(measures),hw1);
end

% Delete Waitbairs
hw1.delete();
hw2.delete();
hw3.delete();

% save SIM_13_Nodes_den_0.1.mat

