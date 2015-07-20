%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: shiftXBs                                                  %
%   Date Started: 9/18/2008                                             %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function accepts the struct params (containing all
%   model parameters), a steady-state value Ca (Ca_max), and a new value 
%   of XB duty cycle.  The AeRU matrix is updated using the new d_cycle, and
%   a new steady state vector of state vars (x0) is computed.  This
%   simulates the rapid increase in detachment (g) that ostensibly occurs
%   during the slack/restretch maneuvre.  It is assumed that the thin
%   filament comes into rapid equilibrium with the new population of XBs
%   (now cycling at a much lower duty cycle).  This new population attaches
%   to the thin filament in a different spatial pattern that influences the
%   dynamics of force redevelopment.
%
%   4/14/09 - Updated function to work correctly with no-XB binding eRU
%   model.
%
%   3/13/10 = Updated function to use the 3-state XB model
%***********************************************************************%

function x0 = shiftXBs(params,      ... % Parameter struct containing all model parameters
                       Ca_i,        ... % Free Ca concentration
                       mult_ktr)       % multiplier of isometric g_app and g_xb to use to simulate ktr
                   

%------------------------%
% Update AeRU with new g %
%------------------------%

bparams_new     = params.bparams;
g_new           = mult_ktr * bparams_new.g;
gxb_new         = mult_ktr * bparams_new.gxb;
bparams_new.g   = g_new;
bparams_new.gxb = gxb_new;

iRUparams_temp = getiRUParams(bparams_new);

params.AiRU    = iRUparams_temp.AiRU;
params.bparams = bparams_new; % Update bparams with version containing g_new (just in case...)

%--------------------------%
% Solve for SS using new g %
%--------------------------%

x0 = solveSingleFullSS(params, Ca_i);

return