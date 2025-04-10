% Author: Eduard Clotet
% Robotics Laboratory
% Universitat de Lleida

function matrix2WordTable(matrixIn, varargin)
%matrix2WordTable(M)
% Copy the content of a Matlab matrix to a table in a Word document using the
% Windows clipboard.
%   Copy matrix into Word Table
%       1) Run the matrix2WordTable function.
%       2) Place the cursor at the position where you want the table to
%            be created inside the Word document.
%       3) Wait for the table to be created and filled automatically.
%
% If values are not being copied properly, increase the CopyPause time.
%           Ex: matrix2WordTable(M, 'CopyPause', 0.7);
%
%**************************************************************************
%Syntax:
%**************************************************************************
%   matrix2WordTable(M)
%   matrix2WordTable(M, Name, Value, ...)
%
% Input Arguments:
%   (Required)
%   M               - The matrix to be copied.
%
% Name-Value Pair Arguments:
%   (Optional)
%   InitPause       - Time (in seconds) given to the user to navigate to 
%                     the word document and select the eqation field.
%                     (5 seconds by default)
%
%   CopyPause       - Time (in seconds) given to MS word to copy a string.
%                     (0.3 seconds by default)
%
%   CreateSpace     - Flag to indicate if the function needs to create the
%                     matrix within the Microsoft Equation field or table
%                     in the Word document.
%                     [0 | 1(default)]
%
%   Format          - Specify the output format of the copied numeric values.
%                     To see the description of this field in detail, check
%                     the formatSpec section in num2str documentation.
%                     
%                     To set the same format to all cells of the matrix:
%                       - formatStr 
%                       *where:
%                           - formatStr is a valid formatSpec string.
%
%                     To specify the format for a certain range of cells:
%                       - {rowsIndx, colsIndx, formatStr}
%                       *where:
%                           - rowsIndx is the index of the rows whose
%                             format will be modified (empty for all rows).
%                           - columnsIndx is the index of the columns whose
%                             format will be modified (empty for all cols).
%                           - formatStr is a valid formatSpec string.
%
%   DecimalMarker   - Single character used to notate the decimal part of 
%                     numbers.
%                     ('.' by defalut)
%
%
%**************************************************************************
%Examples:
%**************************************************************************
%   - Example 1: Create and fill a Word table
%       1) Execute in a Matlab terminal:
%           matrix2WordTable(rand(3,3));
%       2) Navigate to the Word document and place the cursor where the
%          table must be created.
%       3) After 5 seconds (by default) the values will start to get copied
%           to the table.
%
%   - Example 2: Fill the values of an already existing Word Table
%       1) Create a table in the Word document.
%           MS Word menu bar -> insert -> table (select 6x4 size)
%       2) Execute in a Matlab terminal:
%           matrix2WordTable(rand(4,6), 'CreateSpace', 0);
%       3) Navigate to the Word document and select the top-left cell of
%           the table.
%       4) After 5 seconds (by default) the values will start to get copied
%           to the table.
%    
%   - Example 3: Create and fill a Word table with numbers containing four
%	  decimal places:
%       1) Execute in a Matlab terminal:
%           matrix2WordTable(rand(3,3),'Format','%.4f');
%       2) Navigate to the Word document and place the cursor where the
%          table must be created.
%       3) After 5 seconds (by default) the values will start to get copied
%           to the table.

    % Obtain matrix size.
    [numRows, numCols] = size(matrixIn);
    
    % Generate options structure.
    options = parseOptions(numRows, numCols, varargin);
    
    % Intialize a COM Automation server instance with Windows script Host.
    try
        h = actxserver('WScript.Shell');
    catch
        error('Error: Unable to initialize WScript');
    end
    
    % Give the user some time to navigate to the Word document and place
    % the cursor on the equation field.
    disp(['Matlab will begin to copy the values in ',num2str(options.initPause),...
        ' seconds, place the cursor over the MS Equation box!']);
    pause(options.initPause);
    
    % Check if the user wants to create a new matrix
    if options.createSpace
        createTableSpace(h, numCols, options);
    end

    % Create a cell array containing the string that will be pasted to each
    % cell of the matrix
    values = parseMatrix(matrixIn, numRows, numCols, options);
        
    % Use clipboard to copy/paste values to the Word Equation Matrix
    dumpValues(h, values, numCols, options);
    
    %*********************************************************************
    % METHODS
    %*********************************************************************
    function dumpValues(h, values, numCols, options)
        %DUMPVALUES copies the formated strings stored in "values" into the
        %corresponding cells of a Word Equation Matrix through the clipboard.
        %
        % Input parameters:
        % - h: Windows Script Host handle
        % - values: cell array containing the strings to be copied
        % - numCols: number of columns in the table/equation
        % - options: structure containing the parsed options
        
        for vID = 1:length(values)
            clipboard('copy',values{vID});
            h.SendKeys('^{v}');
            pause(options.copyPause);
            h.SendKeys('{RIGHT}');
            pause(options.copyPause);
            if mod(vID,numCols) == 0 && vID ~= length(values)
                if options.createSpace
                    h.SendKeys('{ENTER}');
                else
                    h.SendKeys('{Right}');
                end
                pause(options.copyPause);
            end
        end
    end

    function createTableSpace(h, numCols, opt)
        %CREATETABLESPACE Creates the table where data will be pasted.
        %
        % Input parameters:
        % - h: Windows Script Host handle
        % - numRows: number of rows of the table to be copied
        % - numCols: number of columns of the table to be copied
        % - opt: options structure
        
        str = '+';
        for colID = 1:numCols
            str = strcat(str,{' '},'+');
        end
        
        clipboard('copy',str{1});
        h.SendKeys('^{v}','true');
        pause(opt.copyPause)
        h.SendKeys('{ENTER}','true');
    end
    
    function options = parseOptions(numRows, numCols, params)
        %PARSEOPTIONS initializes the options structure. 
        %
        % Input parameters:
        % - numRows: number of rows of the matrix to be copied
        % - numCols: number of columns of the matrix to be copied
        % - params: user-specified parameters
        % 
        % Output parameters:
        % - options: options strucutre containing;
        %       - formatsStr: cell array of strings containing the format
        %           of each cell.
        %       - copyPause: pause between pastes (in seconds).
        %       - createMatrix: flag that specifies if we need to create
        %           the matrix.
        %       - decimalMarker: marker used to notate the decimal part.

        options.copyPause = 0.3;
        options.createSpace = 1;
        options.decimalMarker = '.';
        options.formatsStr = cell(numRows, numCols);
        options.initPause = 5;
        
        errorFlag = 0;
        
        options.formatsStr = editFormat('%0.3f', numRows,...
            numCols, options.formatsStr);
                    
        for inputID = 1:2:length(params)
            fieldName = params{inputID};
            fieldValue = params{inputID+1};
            
            switch lower(fieldName)
                case 'format'
                    options.formatsStr = editFormat(fieldValue, numRows,...
                        numCols, options.formatsStr);
                case 'initpause'
                    if ~isnumeric(fieldValue) || fieldValue <= 0
                        errorFlag = 1;
                    end
                    options.initPause = fieldValue;
                case 'copypause'
                    if ~isnumeric(fieldValue) || fieldValue <= 0
                        errorFlag = 1;
                    end
                    options.copyPause = fieldValue;
                case 'createspace'
                    if ~isnumeric(fieldValue) || (fieldValue ~= 1 && fieldValue ~= 0)
                        errorFlag = 1;
                    end                    
                    options.createSpace = fieldValue;
                case 'decimalmarker'
                    if ~ischar(fieldValue)
                        errorFlag = 1;
                    end
                    options.decimalMarker = fieldValue;
                otherwise
                    errorFlag = 1;
            end
            
            if errorFlag == 1
                error(['Invalid pair key-value for field ',fieldName]);
            end
        end
        
        %******************************************************************
        % Internal functions
        %******************************************************************
        function formatsStr = editFormat(fieldValue, numRows, numCols, formatsStr)
            if ischar(fieldValue)
                formatsStr(:) = {fieldValue};
            else
                if ~iscell(fieldValue) || length(fieldValue) ~= 3
                    error('Error: Invalid format');
                end
                
                rowIndexs = fieldValue{1};
                colIndexs = fieldValue{2};

                if isempty(colIndexs)
                    colIndexs = 1:numCols;
                end

                if isempty(rowIndexs)
                    rowIndexs = 1:numRows;
                end

                numCols = length(colIndexs);
                numRows = length(rowIndexs);
                for rowID = 1:numRows
                    for colID = 1:numCols
                        formatsStr{rowIndexs(rowID), colIndexs(colID)} = fieldValue{3};
                    end
                end
            end
        end
    end

    function values = parseMatrix(matrixIn, numRows, numCols, options)
        %PARSEMATRIX parses the matrix to a cell array containing the
        %string that will be pasted to each cell of the Word matrix,
        %accounting for format specifications.
        %
        % Input parameters:
        % - matrixIn: user-defined matrix.
        % - numRows: number of rows of the matrix to be copied
        % - numCols: number of columns of the matrix to be copied
        % - opt: options structure
        
        values = cell(1,numRows*numCols);
        linealIndx = 1;
        for rID = 1:numRows
            for cID = 1:numCols
                format = options.formatsStr{rID, cID};
                if isempty(format)
                    values{linealIndx} = num2str(matrixIn(rID,cID));
                else
                    values{linealIndx} = num2str(matrixIn(rID,cID), format);
                end
                
                if options.decimalMarker ~= '.'
                    values{linealIndx} = strrep(values{linealIndx},'.',options.decimalMarker);
                end
                
                linealIndx = linealIndx+1;
            end
        end
    end
end