%********************************************************************%
%   MFM-l_dep
%   Function: calcInitialSL.m
%   Author: Stuart Campbell
%   Date Started: 01/03/2011
%
%   Description: This function estimates the SL that will balance
%   steady-state force between the active and parallel forces (which sum)
%   and the force from the series elastic element.
%
%   NOTE THAT THIS SHOULD NOT BE CALLED IN THE ISOMETRIC CASE (CONSTANT
%   SARCOMERE LENGTH)!!!
%
%********************************************************************%

function SL_init = calcInitialSL(x_initial, ...   The vector of initial conditions for the model
                                 bparams,   ...   The parameter struct
                                 l_tot_initial) % The initial value of the total muscle length
                             
% Use fminsearch

SL0 = 2.0;  % Initial guess

SL_init = fminsearch(@(SL) residFunc(SL, x_initial, bparams, l_tot_initial), SL0);

return


% Helper function to calculate the residual
function resid = residFunc(SL, x_initial, bparams, l_tot_initial)

x_initial(end-2) = SL;  % Replace SL slot with guess

Fp = calcFParallel(bparams.k_p, bparams.tau_p, bparams.SL_slack, SL);
Fs = calcFSeries(bparams.k_s, bparams.tau_s, bparams.l_slack, l_tot_initial - SL);

Fact = calcActiveForce(x_initial, bparams.kxb);

resid = (Fs - (Fp + Fact))^2;

return