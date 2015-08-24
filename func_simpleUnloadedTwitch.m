function rawdata = func_simpleUnloaded_Twitch(varargin)
%-------------------------------------------------------------------------%
% func_simpleUnloadedTwitch.m
% 
% Description: This function generates a simple unloaded contraction
% simulation. It accepts a single argument which should be a valid filename
% referring to a .txt file that contains any parameter changes desired that
% are different from those values contained in the getParams.m function.
%
% The format for the text file is to be tab-delimited: 
% <param_name> <space> <param_value> \n
%-------------------------------------------------------------------------%

% Handle defaults %
fnames = struct('paramfname', 'param_file.txt', ...
                'outputfname', 'output_file.txt');
            
input_arg_names = fieldnames(fnames);
for i = 1:nargin
    fnames.(input_arg_names{i}) = varargin{i};
end



%----------------------------%
% Set Basic Model Parameters %
%----------------------------%

numRUs   = 13;          % Number of RUs in eRU 'rings'

t_end    = 500;         % [ms] - Duration of timecourse simulation

Ca_code  = 2; 
CaT_min  = 0.1;         % [µM] - Minimum Ca concentration
CaT_max  = 1.0;         % [µM] - Maximum Ca concentration


%----------------------------%
%  Set Ca and LS drivers     %
%----------------------------%

Ca_params = {Ca_code, CaT_min, CaT_max};
SL_params  = {1 4 4 0 0}; 

%---------------------------------------%
% Generate argkey/argvals from textfile %
%---------------------------------------%

[argkey, argvals] = generateArgChangesFromTextFile(fnames.paramfname);

% make it an unloaded twitch:
argkey{end+1} = 'k_s';
argvals{end+1} = 0;

%-----------------------------------%
%   Get Auto-generated model data   %
%-----------------------------------%

bparams = getParams(argkey, argvals); % Get basic model params
params  = getFullModelParams(numRUs, bparams);      % Get all auto generated structures


%***************************%
%   Run Twitch Simulation   %
%***************************%

tic
rawdata = twitchResponse(params, t_end, Ca_params, SL_params);
toc


% Plot result
h = figure('Visible', 'off');
plot(rawdata.T, rawdata.SL)
xlabel('Time (ms)')
ylabel('Sarcomere Length (µm)')
print(h, '-dpng', 'output_plot.png')


% Save data to disk
colheads = fieldnames(rawdata); % Get Column headings (fields) of rawdata
colheadmap = true(length(colheads),1);
x_dex      = findCellEntriesContainingTargetStr(colheads,'X');
colheadmap(x_dex) = 0;  % Exclude 'X' from output (it's giant and won't mean much to most people)
colheads   = colheads(colheadmap);

% Loop over rawdata fields and append each column to a cell array of table
% entries
output_cell = cell(length(rawdata.T), length(colheads));    % Preallocate output cell array

for i = 1:length(colheads)
    field = colheads{i};
    output_cell(:,i) = num2cell(rawdata.(field));
       
end

% Add column headings
output_cell = [colheads'; output_cell];

dlmcell(fnames.outputfname, output_cell)
    