%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Script: script_endoCellShorteningMatch.m                            %
%   Date Started: 01/03/2011
%   Author: Stuart Campbell
%                  
%   Description: Simulate a simple unloaded shortening and tweak params to
%   make it match a specific experimental measurement.
%***********************************************************************%

cls

savedir = '/Users/stuart/MATLAB_files/Fellowships2010/';
           

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
SL_slack = 1.66;        % Basically determines the starting or diastolic SL (there will be some active tension at diastolic Ca, which will make the resting value smaller)
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

sim = rawdata.unloaded;
figure(1)
set(gcf, 'Units', 'inches', 'position', [-12.6 3 10.5 3], 'color', 'w')
subplot(141)
hold on
plot(sim.T, sim.SL)

subplot(142)
hold on
plot(sim.T, vnorm(sim.SL))

subplot(143)
hold on
plot(sim.T, sim.X(:,4:5))
title('XB state occupancy')

subplot(144)
hold on
plot(sim.T, sim.X(:,end-1:end))

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

% Calculate quality of the fit
r = calcRsq(t_sl, sl, sim.T, sim.SL);

%-------------------------------%
% Save the Results for Plotting %
%-------------------------------%

sim_t   = sim.T;
sim_SL  = sim.SL;
sim_Ca  = sim.Ca;
data_t  = t_sl;
data_SL = sl;

save([savedir 'endoShorteningMatchData.mat'], 'sim_t', 'sim_SL', 'sim_Ca', 'data_t', 'data_SL')

%%
% Make a parameter table
names = fieldnames(bparams);
numnames = length(names);


data = cell(numnames+1, 3);
data(1,:) = {'Parameter', 'Value', 'Units'};
for i = 1:numnames
    data(i+1,:) = {names{i}, bparams.(names{i}), 'HERE'};
end

fmt_exceptions = cell(numnames + 1, 3);

fmt_exceptions(4+1, 2)       = {'%.3f'};
fmt_exceptions((5:9)+1, 2)   = {'%.3e'};
fmt_exceptions((13:14)+1, 2) = {'%.1e'};

makeTabDelimTable(data, 'endo_fit_params.txt', fmt_exceptions)

% system('python ../common/tab_to_latex.py endo_fit_params.txt /Users/stuart/gen_docs/MSM_Genetic_HF_Review/DataSupplement/endo_fit_params.tex')
