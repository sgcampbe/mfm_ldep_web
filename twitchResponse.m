%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: twitchResponse                                            %
%   Date Started: 8/26/2008                                             %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function handles simulation of a single twitch,
%   including the obtaining of steady-state initial conditions based on
%   minimum initial Ca concentration.
%   See Program Glossary for variable definitions.
%***********************************************************************%

function data = twitchResponse(params,... % Struct containing equation data
                                t_end, ... % [ms] - duration of simulation
                                Ca_params, ... % Parameters specifying driving Ca
                                SL_params)

%---------------------------%
% Set up Initial Conditions %
%---------------------------% 

Ca_min = getCa(0.1, Ca_params);

x0 = solveSingleFullSS(params, Ca_min);

%----------------%
% Integrate ODEs %
%----------------%

data = solveDynamic(x0, params, 0, t_end, Ca_params, SL_params);


return
