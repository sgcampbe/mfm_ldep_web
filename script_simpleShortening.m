%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Script: script_simpleShortening.m                                   %
%   Date Started: 01/03/2011
%   Author: Stuart Campbell
%                  
%   Description: Simulate a simple unloaded shortening twitch to test
%   Ldep machinery.
%***********************************************************************%

cls

%---------------------------%
%   Simulation Parameters   %
%---------------------------%

numRUs   = 13;          % Number of RUs in eRU 'rings'

t_end    = 2000;        % [ms] - Duration of timecourse simulation

Ca_code  = 10;          % Digitized from ENDO rat cell 

CaT_min  = 0.1;         % [µM] - Minimum Ca concentration
CaT_max  = 1.0;         % [µM] - Maximum Ca concentration



%----------------------------------%
%      Length dependent params     %
%----------------------------------%

k_s      = 0.0;         % Setting this to zero causes the cell to shorten against only the internal load (k_p)
SL_slack = 1.66;       % Basically determines the starting or diastolic SL (there will be some active tension at diastolic Ca, which will make the resting value smaller)
k_p      = 0.01;
tau_p    = 5;
kxb      = 10.8;
nu       = 0.003;

% Tweak other params here, if desired:
ldep_key  = {'k_s', 'SL_slack', 'k_p', 'tau_p', 'kxb', 'nu'};
ldep_vals = { k_s,   SL_slack,   k_p,   tau_p,   kxb,   nu };

% Set up L_params
L_code = 1; % Isometric total muscle length - just a formality since series stiffness is zero
L_min  = 4; % Doesn't really matter - just needs to be longer than SL_slack, probably


%----------------------------------%
% Parameterization of RLCP effects %
%----------------------------------%

rlcpkey  = {'p_f', 'p_hf', 'p_hb', 'p_g', 'p_x0'};
rlcpvals = { 0.6,   -0.4,   0.4,    0.3,    0.3 }; 

% Add other mutual parameter values
basekey  = [rlcpkey  {'Temp', 'fapp', 'gxb', 'xbmodsp', 'kon', 'gammaB', 'RLCP', 'ku'}];
basevals = [rlcpvals {273+25,  0.20,     0.225,   1.0,     350,   300,      0.3,   0.35}];

%-------------------------------%
% Condition-specific Parameters %
%-------------------------------%

% Unloaded Shortening
argkey.unloaded   = [basekey  ldep_key]; 
argvals.unloaded  = [basevals ldep_vals];     


names = fieldnames(argkey);

%----------------------%
% Loop over conditions %
%----------------------%

for i = 1:length(names)
    
    name = names{i};
    
    %----------------------------%
    %  Set Ca Transient driver   %
    %----------------------------%

    Ca_params = {Ca_code(i) CaT_min CaT_max};
    L_params  = {L_code L_min L_min 0 0};

    %-----------------------------------%
    %   Get Auto-generated model data   %
    %-----------------------------------%

    bparams = getParams(argkey.(name), argvals.(name)); % Get basic model params
    params  = getFullModelParams(numRUs, bparams);      % Get all auto generated structures


    %***************************%
    %   Run Twitch Simulation   %
    %***************************%

    tic
    rawdata.(name) = twitchResponse(params, t_end, Ca_params, L_params);
    toc
    
    
    %--------------------------%
    %   Data Post-processing   %
    %--------------------------%
    % rawdata.(name).force = calcForce(rawdata.(name).X, params, rawdata.(name).T, SL_params);
    
    % % Calculate properties of twitches and model fit
    % props.(name)  = transientProps(rawdata.(name).T, rawdata.(name).force);

end


%------------------------%
%   Plot Final Results   %
%------------------------%

data = rawdata.unloaded;
figure(1)
set(gcf, 'Units', 'inches', 'position', [-12.6 3 10.5 3], 'color', 'w')
subplot(141)
hold on
plot(data.T, data.SL)

subplot(142)
hold on
plot(data.T, vnorm(data.SL))

subplot(143)
hold on
plot(data.T, data.X(:,4:5))
title('XB state occupancy')

subplot(144)
hold on
plot(data.T, data.X(:,end-1:end))

% Load some data for fitting purposes
data_dir = '/Users/stuart/data/PreliminaryMyocyteData/';
load([data_dir '101222_rat/processed_data.mat']);
t_sl = 1000 * avg_data.endo(1).time;    % Convert to ms
sl = avg_data.endo(1).SL;
subplot(141)
hold on
plot(t_sl, sl, 'r')

subplot(142)
hold on
plot(t_sl, vnorm(sl), 'r')