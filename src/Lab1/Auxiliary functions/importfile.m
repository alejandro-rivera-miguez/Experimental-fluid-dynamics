%% =======================================================================
% FUNCTION: importfile
% DESCRIPTION: Automatically generated function to read raw calibration 
%              voltage data from the external force balance.
%
% AUTHORS:  Alejandro Rivera Míguez
%           Alberto Rivero García
%           Mikel Segovia Díaz
% =======================================================================
function Taraturabilanciabolla = importfile(filename, dataLines)

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
    opts.VariableNames = ["NI9237ch0VVch1VVch2VVNI9205Ch0Vch1Vch2V", "VarName2", "VarName3", "Var4", "Var5", "Var6"];
    opts.SelectedVariableNames = ["NI9237ch0VVch1VVch2VVNI9205Ch0Vch1Vch2V", "VarName2", "VarName3"];
    opts.VariableTypes = ["double", "double", "double", "string", "string", "string"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Specify variable properties
    opts = setvaropts(opts, ["Var4", "Var5", "Var6"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["Var4", "Var5", "Var6"], "EmptyFieldRule", "auto");
    
    % Import the data
    Taraturabilanciabolla = readtable(filename, opts);
end