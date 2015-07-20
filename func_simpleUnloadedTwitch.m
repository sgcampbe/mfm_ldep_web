function rawdata = func_simpleUnloaded_Twitch(paramfname)
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

[argkey, argvals] = generateArgChangesFromTextFile(paramfname);

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
plot(rawdata.T, rawdata.SL)