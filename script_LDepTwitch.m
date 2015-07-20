%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Script: script_LDepTest.m (from script_WTDM2hz4hzTwitches.m         %
%   Date Started: 7/12/2010 (Jared) 
%   Author: Stuart Campbell and Jared Tangney                           %
%                                                                       %
%   Description: This script utilizes a dynamic SL to predict force 
%   outputs.
% **********************************************************************%

clear all

%---------------------------%
%   Simulation Parameters   %
%---------------------------%

numRUs   = 13;          % Number of RUs in eRU 'rings'

%SL       = 2.0;         % Note bogus SL - not yet implemented!

t_end    = 500;         % [ms] - Duration of timecourse simulation

Ca_code  = [2]; 
%Ca_code  = [7 8 7 8 8 8 8 8 7];   % Type of Ca input, corresponds to measured CaT at 4 Hz (code 8) and 2 Hz (code 7)
CaT_min  = 0.1;         % IGNORED! [µM] - Minimum Ca concentration
CaT_max  = 1.0;         % IGNORED! [µM] - Maximum Ca concentration



%----------------------------------%
%      Length dependent params     %
%----------------------------------%



%----------------------------------%
% Parameterization of RLCP effects %
%----------------------------------%

rlcpkey  = {'p_f', 'p_hf', 'p_hb', 'p_g', 'p_x0'};
rlcpvals = { 0.6,   -0.4,   0.4,    0.3,    0.3 }; 

% Add other mutual parameter values
rlcpkey  = [rlcpkey  {'Temp', 'fapp', 'gxb', 'xbmodsp', 'kon', 'gammaB'}];
rlcpvals = [rlcpvals {273+25,  0.20,     2.0,   1.0,     350,   300}];

%-------------------------------%
% Condition-specific Parameters %
%-------------------------------%

% WT at 2 Hz 
argkey.WT2   = [rlcpkey  {'RLCP', 'ku'}]; 
argvals.WT2  = [rlcpvals { 0.3,   0.45}];     


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
    SL_params = {1 2.0 2.1 150 150};                  % SL_code SLmax SLmin iso_t1 strDur

    %-----------------------------------%
    %   Get Auto-generated model data   %
    %-----------------------------------%

    bparams = getParams(argkey.(name), argvals.(name)); % Get basic model params
    params  = getFullModelParams(numRUs, bparams);      % Get all auto generated structures


    %***************************%
    %   Run Twitch Simulation   %
    %***************************%

    tic
    rawdata.(name) = twitchResponse(params, t_end, Ca_params, SL_params);
    toc
    
    
    %--------------------------%
    %   Data Post-processing   %
    %--------------------------%
    rawdata.(name).force = calcForce(rawdata.(name).X, params, rawdata.(name).T, SL_params);
    
    % Calculate properties of twitches and model fit
    props.(name)  = transientProps(rawdata.(name).T, rawdata.(name).force);

end


%------------------------%
%   Plot Final Results   %
%------------------------%

%subplot(211)
%for i = [1 2 5 6 7]
%    name = names{i};
%    plot(rawdata.(name).T, rawdata.(name).force)
%    hold on
%end

%subplot(212)
%for i = [3 4 8]
%    name = names{i};
%    plot(rawdata.(name).T, rawdata.(name).force)
%    hold on
%end

figure(1)
plot(rawdata.WT2.T, rawdata.WT2.force, 'c')
hold on
%plot(rawdata.WT2.T, SL, 'k')

%plot(rawdata.DM2_notnip.T, vnorm(rawdata.DM2_notnip.force), 'r')



%------------------------%
%   Save Final Results   %
%------------------------%

simdata = rawdata;  % rename...
save SimData/WTDM2hz4hzTwitches.mat simdata props 