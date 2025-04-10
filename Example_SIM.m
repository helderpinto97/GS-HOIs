%% Example Code ---  Greedy Approach to Quantify Redundancy/Synergy
clear all; close all; clc;

% Addpath for functions
addpath('functions\');

% Measures that the framework support
measures={'iir','rsi','oir'};

% Fix Seed for Reproducibility
rng('default');

% Parameters
ii=3; % Target Process
model_order_p=1; % Order of the Model
q=22; % Auto correlation truncation lag
numsurr=100; % number of surrogates
numsim=100; % number of realizations for each coupling parameter
alpha= 0.05; % Significance level for the surrogate test

% Coupling parameter (X1 and X3)
c=0:0.1:1;

hw1= waitbar(0,'Simulation Loop...');
hw2= waitbar(0,'Coupling Parameter...');
hw3=waitbar(0,'Measure...');
for meas=1:length(measures)
    for ic=1:length(c)
        for isim=1:numsim
            % Noise Vector
            x=randn(1000,4);

            % Equations X_3 and X_4
            for t=2:1000
                x(t,4)=0.9*x(t-1,3)+0.1*x(t,4);
                x(t,3)=c(ic)*x(t-1,1)+0.5*x(t-1,2)+0.1*x(t,3);
            end
            % Series on row
            x=normalize(x);
            Y=x';

            % Greedy Approach 
            ret=greedy_red_syn(Y,model_order_p,ii,q,numsurr,alpha,measures{meas});

            if meas==1
                Redundancy_est(isim,ic,meas)=ret.IIRM(end);
                Redundancy_cM{isim,ic,meas}=ret.IIRM_cond_vec;
                Redundancy_in_triplet{isim,ic,meas}=ret.IIRM_in_triplet;
                Synergy_in_triplet{isim,ic,meas}=ret.IIRm_in_triplet;
                Synergy_est(isim,ic,meas)=ret.IIRm(end);
                Synergy_cm{isim,ic,meas}=ret.IIRm_cond_vec;
            elseif meas==2
                Redundancy_est(isim,ic,meas)=ret.RSIM(end);
                Redundancy_cM{isim,ic,meas}=ret.RSIM_cond_vec;
                Redundancy_in_triplet{isim,ic,meas}=ret.RSIM_in_triplet;
                Synergy_in_triplet{isim,ic,meas}=ret.RSIm_in_triplet;
                Synergy_est(isim,ic,meas)=ret.RSIm(end);
                Synergy_cm{isim,ic,meas}=ret.RSIm_cond_vec;
            else
                Redundancy_est(isim,ic,meas)=ret.OIRM(end);
                Redundancy_cM{isim,ic,meas}=ret.OIRM_cond_vec;
                Redundancy_in_triplet{isim,ic,meas}=ret.OIRM_in_triplet;
                Synergy_in_triplet{isim,ic,meas}=ret.OIRm_in_triplet;
                Synergy_est(isim,ic,meas)=ret.OIRm(end);
                Synergy_cm{isim,ic,meas}=ret.OIRm_cond_vec;
            end
            waitbar(isim/numsim,hw1);
        end
        waitbar(ic/numel(c),hw2);
    end
    waitbar(meas/numel(measures),hw3);
end

% Delete waitbars
hw1.delete();
hw2.delete();
hw3.delete();


%% Plot Results
measures_s_r={'Synergy_est','Redundancy_est'};
col={'b','r'};
for meas=1:length(measures)
    fig= figure('WindowState','maximized','Name',measures{meas});
    for j=1:length(measures_s_r)
        s=shadedErrorBar(c,eval([measures_s_r{j} '(:,:,meas)']),{@median,@(x) dist_qtl(x)},'lineprops',col{j},'transparent',true,'patchSaturation',0.075);
        s.mainLine.LineWidth = 1.5; hold on;
        % plot(c,eval(measures_th{j}),'Color', 'k', 'LineWidth',1.5,'DisplayName','Theoretical Value');
        axis 'square';
        box on;
        ylabel('$\Phi_{IIR}$','Interpreter','latex','FontName','Helvetica');
        xlabel('$c$','Interpreter','latex','FontName','Helvetica');
        legend({'$\Phi_{IIR}$ - Synergy','$\Phi_{IIR}$ - Redundancy'},'Box', 'off','Interpreter','latex','Location','northwest');
        ax=gca;
        ax.FontSize=18;
        % title(measures_title{j},'Interpreter','latex','FontWeight','bold');
    end
end
