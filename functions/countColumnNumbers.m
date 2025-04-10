function counts = countColumnNumbers(matrix)
    % This function counts the appearance of numbers in each column of a matrix.
    % Input:
    %   matrix - A numeric matrix
    % Output:
    %   counts - A cell array with tables for each column showing unique numbers and their counts.

    % Get the number of columns
    [~, numCols] = size(matrix);
    
    % Initialize the cell array to hold the counts
    counts = cell(1, numCols);
    
    % Loop through each column
    for col = 1:numCols
        % Get the unique numbers and their counts for the column
        [uniqueNumbers, ~, idx] = unique(matrix(:, col));
        countsInCol = accumarray(idx, 1);
        
        % Store the results in a table for better readability
        counts{col} = table(uniqueNumbers, countsInCol, ...
            'VariableNames', {'Number', 'Count'});
    end
end

