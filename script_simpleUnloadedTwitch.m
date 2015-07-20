%-------------------------------------------------------------------------%
% script_simpleUnloadedTwitch.m
% 
% Description: This script illustrates a simple unloaded contraction
% simulation.
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

%-------------------------------------------%
% Some reasonable parameters for 25 C mouse %
%-------------------------------------------%

argkey  = {'p_f', 'Temp', 'fapp', 'gxb', 'xbmodsp', 'kon', 'gammaB'};
argvals = {0.6    273+25,  0.20,     2.0,   1.0,     350,   300};

%-------------------------------%
% Condition-specific Parameters %
%-------------------------------%

argkey   = [argkey  {'RLCP', 'ku', 'k_s'}];    %Note that k_s = 0 makes unloaded shortening 
argvals  = [argvals { 0.3,   0.45,   0}];    

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