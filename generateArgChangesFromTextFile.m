function [argkey, argvals] = generateArgChangesFromTextFile(paramfname)

%-------------------------------------------------------------------------%
% This function takes a space-delimited text file where each line contains a
% paramter name <space> parameter value pair and parses it into the cell
% arrays argkey and argvals, ready to be fed into the getParams function.
%-------------------------------------------------------------------------%


%--------------------%
% Read in param file %
%--------------------%

paramstr = file2str(paramfname);
lines    = split(paramstr, '\n');


%-------------------------------%
% Extract values from each line %
%-------------------------------%

numparams = length(lines);

% Initialize variables
argkey   = cell(1,numparams);
argvals  = cell(1,numparams); 

for i = 1:numparams
    parts = split(lines{i});    % Split by whitespace
    argkey{i} = parts{1};
    argvals{i} = str2double(parts{2});
end
