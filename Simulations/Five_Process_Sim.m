clear all; close all; clc;
addpath("../functions");

%% Parameters
M=5; %%% Number of Processes
N=5000; %%% Length of Realization
model_order_p=2; % Model Order
tau=ones(1,M); % Time Delay 
q=22; % Auto correlation truncation lag
nsurr=100; % number surrogates
nsim=100; % number of simulations
alpha=0.05;

% X and Y labels for HeatMaps
labels_plot={'$X_1$','$X_2$','$X_3$','$X_4$','$X_5$'};

% Measures that the framework support
measures={'iir','rsi','oir'};


% Target Process
ii=1;

% Coupling 
c=0:0.1:1;

% Waitbars
hw1 = waitbar(0,'Coupling Parameter...');
hw2= waitbar(0,'Simulation Loop...');
hw3= waitbar(0,'Measure Loop...');

for meas=1:length(measures)
    for ic=1:length(c)

        % Poles of Nodes --- Put the initial and central node with HF oscillations
        % High Frequency Nodes
        par.poles{1}=([0.85 0.1]);
        par.poles{2}=([0.85 0.1]);
        par.poles{3}=([0.85 0.1]);
        par.poles{4}=([0.85 0.35]);

        % Noise poles
        par.poles{5}=([0.7 0.2]);
        % par.poles{6}=([]);

        % Coupling
        par.coup=[2 1 1 1-c(ic);
            3 1 2 c(ic);
            4 2 2 1;
            4 3 2 1];

        % Diagonal of Covariance Matrix
        par.Su=ones(1,M); %variance of innovation processes

        % VAR Parameters
        [Am,Su]=theoreticalVAR(M,par); %% VAR parameters

        % Transpose Parameter Matrix
        Am=Am';

        for isim=1:nsim
            U = mvnrnd(zeros(1,M),Su,N);
            Y=var_filter(Am,U);

            % Normalize the Data
            Y=normalize(Y);

            % Series on row
            Y=Y';
            % Greedy Approach
            ret=greedy_red_syn_est(Y,model_order_p,ii,q,nsurr,alpha,measures{meas});

            if strcmp(measures{meas},'iir')
                Redundancy_est(isim,ic,meas)=ret.IIRM(end);
                Redundancy_cM{isim,ic,meas}=ret.IIRM_cond_vec;
                Redundancy_in_triplet{isim,ic,meas}=ret.IIRM_in_triplet;
                Synergy_in_triplet{isim,ic,meas}=ret.IIRm_in_triplet;
                Synergy_est(isim,ic,meas)=ret.IIRm(end);
                Synergy_cm{isim,ic,meas}=ret.IIRm_cond_vec;
            elseif strcmp(measures{meas},'rsi')
                Redundancy_est(isim,ic,meas)=ret.RSIM(end);
                Redundancy_cM{isim,ic,meas}=ret.RSIM_Final_Vec;
                Redundancy_in_triplet{isim,ic,meas}=ret.RSIM_in_triplet;
                Synergy_in_triplet{isim,ic,meas}=ret.RSIm_in_triplet;
                Synergy_est(isim,ic,meas)=ret.RSIm(end);
                Synergy_cm{isim,ic,meas}=ret.RSIm_Final_Vec;
            else
                Redundancy_est(isim,ic,meas)=ret.OIRM(end);
                Redundancy_cM{isim,ic,meas}=ret.OIRM_Final_Vec;
                Redundancy_in_triplet{isim,ic,meas}=ret.OIRM_in_triplet;
                Synergy_in_triplet{isim,ic,meas}=ret.OIRm_in_triplet;
                Synergy_est(isim,ic,meas)=ret.OIRm(end);
                Synergy_cm{isim,ic,meas}=ret.OIRm_Final_Vec;
            end
            waitbar(isim/nsim,hw2);
        end
        waitbar(ic/numel(c),hw1);
    end
    waitbar(meas/numel(measures),hw3);
end

hw1.delete();
hw2.delete();
hw3.delete();



%% Plot Results
measures_s_r={'Synergy_est','Redundancy_est'};
ylabels_title={'$\Phi_{IIR}$ [nats]','$\Phi_{RSIR}$ [nats]','$\Phi_{\Delta \mathrm{OIR}}$ [nats]'};
col={'b','r'};
for meas=1:length(measures)
    fig= figure('WindowState','maximized','Name',measures{meas});
    for j=1:length(measures_s_r)
        s=shadedErrorBar(c,eval([measures_s_r{j} '(:,:,meas)']),{@median,@(x) dist_qtl(x)},'lineprops',col{j},'transparent',true,'patchSaturation',0.075);
        s.mainLine.LineWidth = 1.5; hold on;
        % plot(c,eval(measures_th{j}),'Color', 'k', 'LineWidth',1.5,'DisplayName','Theoretical Value');
        axis 'square';
        box on;
        ylabel(ylabels_title{meas},'Interpreter','latex','FontName','Helvetica');
        xlabel('$c$','Interpreter','latex','FontName','Helvetica');
        legend({'Synergy','Redundancy'},'Box', 'off','Interpreter','latex','Location','northwest');
        ax=gca;
        ax.FontSize=18;
        % title(measures_title{j},'Interpreter','latex','FontWeight','bold');
    end
end
