function esa = import_ESAfile(filename, dataLines)
%IMPORTFILE1 Import data from a text file
%  SWEEPESA2024 = IMPORTFILE1(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  SWEEPESA2024 = IMPORTFILE1(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  SWEEPESA2024 = importfile1("/Users/lykhoo/Library/CloudStorage/OneDrive-PrincetonUniversity/SWAPI_FM_CAL_database/FM_Pre_Calibration/Raw_Data/2024.03.27/passband_3keVHe_IAn90/SWEEP_ESA_2024.03.27_74.csv", [3, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 29-Mar-2024 20:27:14

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [3, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 11);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["TimeStamp", "SupplyName", "CommandSent", "Action", "Voltage", "Current", "HVStatus", "Errors", "PRMRate", "SECRate", "COINRate"];
opts.VariableTypes = ["string", "categorical", "double", "categorical", "double", "double", "categorical", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "TimeStamp", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["TimeStamp", "SupplyName", "Action", "HVStatus"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["CommandSent", "Errors"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["CommandSent", "Errors"], "ThousandsSeparator", ",");

% Import the data
esa = readtable(filename, opts);

end