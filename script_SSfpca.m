%************************************************************************%
% Script: script_ssFpCa.m
% Author: Stuart Campbell
% Date Started: 2/27/10
% Description: This script generates steady-state force-pCa curves at
% a given sarcomere length.
%************************************************************************%

clear

%---------------------------%
%   Simulation Parameters   %
%---------------------------%

numRUs  = 3;           % Number of RUs in eRU model
Ca_min  = pCa2uM(8);    % [µM] - Minimum Ca concentration
Ca_max  = pCa2uM(4);    % [µM] - Maximum Ca concentration
numPts  = 20;           % Number of SS points to simulate

SL      = 2.1;          % Sarcomere length

argkey  = {'gammaB'};
argvals = {30};

%-----------------------------------%
%   Get Auto-generated model data   %
%-----------------------------------%

bparams = getParams(argkey, argvals);                   % Get basic model params
params  = getFullModelParams(numRUs, bparams);  % Get all auto generated structures

%--------------------%
%   Run Simulation   %
%--------------------%

tic
rawdata  = runSSFull(params, Ca_min, Ca_max, numPts, SL);
force_ss = calcForceSS(rawdata.X, params, SL);
toc

%--------------------------%
%   Data Post-processing   %
%--------------------------%

Ca_range = rawdata.Ca_range;

%------------------------%
%   Plot Final Results   %
%------------------------%

fig = figure(1);
hold on
set(fig, ...
    'Units', 'Inches', ...    
    'Position', [7 16 5 5], ...
    'Color', 'w');
for i = 1:length(SL)
    pCaPlot(Ca_range, force_ss(:,i), gca, randColor, '-', 1.25)
    hold on
end