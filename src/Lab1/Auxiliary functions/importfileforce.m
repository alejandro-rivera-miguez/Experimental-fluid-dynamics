%% =======================================================================
% FUNCTION: importfileforce
% DESCRIPTION: Automatically generated function to read raw aerodynamic 
%              voltage data recorded during the wind tunnel testing of 
%              the airfoil at varying angles of attack.
%
% AUTHORS:  Alejandro Rivera Míguez
%           Alberto Rivero García 
%           Mikel Segovia Díaz
% =======================================================================
function gruppo5 = importfileforce(filename, dataLines)

    % If dataLines is not specified, define defaults
    if nargin < 2
        dataLines = [2, Inf];
    end
    
    %% Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 6);
    
    % Specify range and delimiter
    opts.DataLines = dataLines;
    opts.Delimiter = "\t";
    
    % Specify column names and types
    opts.VariableNames = ["Fx", "Fy", "M", "PsmP2", "PtmP2", "P1mP2"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Import the data
    gruppo5 = readtable(filename, opts);
end