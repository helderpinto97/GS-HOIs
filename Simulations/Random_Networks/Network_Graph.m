clear all; close all; clc;

structs={'Model_Struct_13_Nodes_den_0.05','Model_Struct_13_Nodes_den_0.1','Model_Struct_13_Nodes_den_0.25','Model_Struct_13_Nodes_den_0.3'};

for ii=1:length(structs)
    load([structs{ii} '.mat']);

    % Adj Matrix
    A = abs(Data.Model');

    % Create a directed graph from the adjacency matrix
    G = digraph(A);

    % Define node labels
    n = size(A,1); % Define the number of elements you need
    nodesLabels = arrayfun(@(i) sprintf('X_{%d}', i), 1:n, 'UniformOutput', false);

    % Automatically calculate node positions using the 'force' layout
    fig = figure('WindowState','maximized');
    h = plot(G, 'Layout', 'force', 'Marker', 'o', 'MarkerSize', 12, 'NodeColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5, 'ArrowSize', 15);

    % Set edge widths based on adjacency matrix values
    edgeWeights = G.Edges.Weight; % Extract edge weights from the graph
    scaledEdgeWidths = 1 + 4 * (edgeWeights - min(edgeWeights)) / (max(edgeWeights) - min(edgeWeights)); % Scale to range [1,5]
    h.LineWidth = scaledEdgeWidths;

    % Ensure that edges are solid black
    h.EdgeColor = [0, 0, 0];  % Set edge color explicitly to black
    h.EdgeAlpha = 1;  % Ensure the edges are fully opaque

    % Draw perfect circles around node labels and place labels inside
    hold on;
    for i = 1:length(nodesLabels)
        % Get the position of each node from the plot object
        x_pos = h.XData(i);
        y_pos = h.YData(i);

        % Draw circles around node labels and place labels inside
        rectangle('Position', [x_pos-0.2, y_pos-0.2, 0.4, 0.4], 'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceColor', 'w');
        text(x_pos, y_pos, nodesLabels{i}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
    end
    hold off;

    % Set the axes to be equal to ensure perfect circles
    axis equal;
    ax = gca;
    ax.LineWidth = 2;
    ax.FontSize = 20;


    % % Adj Matrix
    % A=abs(Data.Model');
    % 
    % % Create a directed graph from the adjacency matrix
    % G = digraph(A);
    % 
    % % Define node labels
    % n = size(A,1); % Define the number of elements you need
    % nodesLabels = arrayfun(@(i) sprintf('X_{%d}', i), 1:n, 'UniformOutput', false);
    % 
    % 
    % % Automatically calculate node positions using the 'force' layout
    % fig=figure('WindowState','maximized');
    % h = plot(G, 'Layout', 'force', 'Marker', 'o', 'MarkerSize', 12, 'NodeColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5, 'ArrowSize', 10);
    % 
    % % Ensure that edges are solid black
    % h.EdgeColor = [0, 0, 0];  % Set edge color explicitly to black
    % % Remove any alpha transparency (if it exists) and set edge color explicitly
    % h.EdgeAlpha = 1;     % Ensure the edges are fully opaque
    % 
    % % Draw perfect circles around node labels and place labels inside
    % hold on;
    % for i = 1:length(nodesLabels)
    %     % Get the position of each node from the plot object
    %     x_pos = h.XData(i);
    %     y_pos = h.YData(i);
    % 
    %     % Draw circles around node labels and place labels inside
    %     rectangle('Position', [x_pos-0.2, y_pos-0.2, 0.4, 0.4], 'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceColor', 'w');
    %     text(x_pos, y_pos, nodesLabels{i}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
    % end
    % hold off;
    % 
    % % Set the axes to be equal to ensure perfect circles
    % axis equal;
    % ax=gca;
    % ax.LineWidth=2;
    % ax.FontSize=20;
    exportgraphics(fig,[structs{ii} '.eps'], 'ContentType','vector','Resolution',600,'BackgroundColor','none');
end