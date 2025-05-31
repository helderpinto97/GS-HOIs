clear all; close all; clc;

% Addpath 
addpath('../../functions');

% Load results
load('SIM_13_Nodes_den_0.1.mat');
load('Theoretical_Values_13_Nodes_den_0.1.mat');

% Measures 
hois={'Redundancy_est','Synergy_est'};
hois_th={'Redundancy_th','Synergy_th'};
measures={'iir','rsi','oir'};
measures_title={'$\Phi_{IIR}$','$\Phi_{RSIR}$','$\Phi_{\Delta OIR}$'};

N = 13;
X_labels = arrayfun(@(i) sprintf('$X_{%d}$', i), 1:N, 'UniformOutput', false);

fig=figure('WindowState','maximized');
positions_bar=[1 2 3;
    4 5 6];
for i=1:length(hois)
    for j=1:length(measures)
        subplot_tight(2,3,positions_bar(i,j),[0.11 0.05]); hold on;
        axis('square');
        boxplot(eval([hois{i} '(:,:,j)']),'Labels',X_labels,'Colors','k');
        h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
        h=findobj('Marker','+'); set(h, 'Marker','o'); set(h,'MarkerFaceColor','r');
        set(findobj(gca,'type','line'),'linew',2);
        plot(eval([hois_th{i} '(j,:)']),'o','Color',[0.4660 0.6740 0.15],'MarkerFaceColor',[0.4660 0.6740 0.15],'MarkerEdgeColor',[0.4660 0.6740 0.15],'MarkerSize',5,'DisplayName','Theoretical Value');
        if i==1
            title(measures_title{j},'Interpreter','latex','FontWeight','bold');
        end
        % sgtitle(pairs_labels{i},'FontWeight','Bold','Interpreter','latex');

        if i==1 && j==1
            ylabel('Redudancy','FontWeight','bold');
        elseif i==2 && j==1
            ylabel('Synergy','FontWeight','bold');
        end

        if i==1
            set(gca,'xtick',[]);
            set(gca,'xticklabel',[]);
            % ylim([-0.1 0.35]);
            if j==1
                % ylim([0.1 0.4]);
            elseif j==2
                % ylim([1.3 1.9]);
            else
                % ylim([0.3 0.5]);
            end
        else
            ylim([-0.15 0.1]);
        end
        xtickangle(-45);
        if j==3 && i==1
            legend('Box','off','Location','best');
        end
        ax=gca;
        ax.FontSize=18;
        ax.LineWidth=2;
        ax.XAxis.TickLabelInterpreter = 'latex';
        
    end
end

% exportgraphics(fig,'Random_Network_Results_13_Nodes_den_0.1_th.png','Resolution',600,'BackgroundColor','none','ContentType','vector');