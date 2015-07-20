function [err comp_vals] = runTwitchSim(opt_params,    ... Vector of parameters from optimization function
                                         opt_names,     ... Cell array of parameter names corresponding to entries in opt_params
                                         varargin)        % Optional changes to task parameters such as the objective function handle, length of the simulation, etc
                                                          % Will be read in using parse_pv_pairs.
                                                       
%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: runTwitchSim                                              %
%   Date Started: 12/29/2011                                            %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This is intended to be a fairly generic function for
%   evaluating twitches based on parameters passed in, such as those
%   provided by an opmitization routine.  The optimization parameters are
%   passed in via the vector opt_params.  opt_names will match these
%   parameters with those in getParams so that they can overwrite defaults.
%   The length of the twitch simulation is determined by task_params.t_end.
%   The objective function, which is responsible for returning the error
%   (err), is flexibly specified by the function handle in
%   task_params.obj_func_handle.  The objective function must return two
%   arguments: err (the fit error) and comp_vals (a variable, flexibly
%   defined, which somehow holds the simulation properties and data
%   properties, with labels if needed, to show the numbers that went into
%   computing the error.  comp_vals will make it much easier and less
%   error-prone to visualize the fit.  task_params must also hold a field
%   pointing to the .mat file containing data to which the simulation will
%   be compared.  The procedure for making the comparison is determined
%   entirely by the contents of the objective function.  Arguments to the
%   objective function are: t (the simulation time points), SL (simulated
%   sarcomere length), force (simulated force), and task_params.  The 
%   objective function is also an ideal place to perform plotting commands 
%   specific to the fit, so that the code of runTwitchSim is generic.
%   task_params.param_file_name must point to a valid parameter file.  The
%   parameter file defines ranges of parameter values to turn the scaled
%   opt_params values into absolute values.  These are extracted according
%   to opt_names.  If the parameter file contains entries not found in
%   opt_names, then these are considered default values - in other words,
%   they are specified as different than the defaults in getParams, but
%   they are not actively fitted.  
%***********************************************************************%


%----------------------------%
% Set up default task_params %
%----------------------------%

task_params.numRUs   = 13;          % Number of RUs in eRU 'rings'

task_params.t_end    = 2000;        % [ms] - Duration of timecourse simulation

task_params.CaT_min  = 'none';      % [µM] - Minimum Ca concentration
task_params.CaT_max  = 'none';      % [µM] - Maximum Ca concentration

task_params.L_code   = 1;           % Isometric total muscle length - just a formality since series stiffness is zero
task_params.L_min    = 4;           % Doesn't really matter - just needs to be longer than SL_slack, probably

task_params.CaT_file_name   = 'fit_test_data/testCaT.mat';       % String specifying file from which to load the input Ca transient
task_params.data_file_name  = 'fit_test_data/testSHprops.mat';   % String specifying file from which to load data for the fit calculation
task_params.param_file_name = 'fit_test_data/testParamFile.txt'; % String specifying file from which to load parameter ranges and default parameters

task_params.obj_func_handle = @objFunc_test;                     % Handle to the objective function used to calculate error between fit and data
task_params.interim_results_file_handle = [];                    % Handle to file for writing intermediate results - only used if supplied with a valid handle
task_params.best_fit_fname  = 'temp.mat';                        % String specifying file name/location for writing best-so-far info.  Used in plotting in objective function

%---------------------------%
% Unpack passed task_params &
%---------------------------%

task_params = parse_pv_pairs(task_params,varargin);


%-----------------------------------------%
% Store Current opt_params in task_params %
%-----------------------------------------%

task_params.opt_params = opt_params;    % Used in objective function


%----------------------------------------------------------------------%
% Transform opt_params into argvals ready for getParams call; retrieve %
% other defaults from file                                             %
%----------------------------------------------------------------------%

[abs_argvals abs_argkey] = optParamsToAbsParams(opt_params, opt_names, task_params.param_file_name);


%----------------------------%
%  Set Ca Transient driver   %
%----------------------------%

Ca_params = {0, task_params.CaT_min, task_params.CaT_max, task_params.CaT_file_name};
L_params  = {task_params.L_code task_params.L_min task_params.L_min 0 0};


%-----------------------------------%
%   Get Auto-generated model data   %
%-----------------------------------%

bparams = getParams(abs_argkey, abs_argvals, true);         % Get basic model params
params  = getFullModelParams(task_params.numRUs, bparams);  % Get all auto generated structures


%***************************%
%   Run Twitch Simulation   %
%***************************%

tic
rawdata = twitchResponse(params, task_params.t_end, Ca_params, L_params);
toc


%--------------------------%
%   Data Post-processing   %
%--------------------------%
force = calcForce(rawdata.X, params, rawdata.T, L_params);
t     = rawdata.T;
SL    = rawdata.SL;


%-------------------------%
% Call Objective Function %
%-------------------------%

[err comp_vals] = task_params.obj_func_handle(t, SL, force, task_params);


return



