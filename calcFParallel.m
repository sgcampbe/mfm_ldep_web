%************************************************************************%
% Myofilament Model - Length Dependent
% Function: calcFParallel.m
% Author: Stuart Campbell
% Start date: 12/29/10
% Description: Calculate the force due to parallel elastance after the
% formulation of Rice (F_titin in his work)
%************************************************************************%

function F_parallel = calcFParallel(k_p, tau_p, SL_slack, SL_curr)

if SL_curr >= SL_slack
    F_parallel = k_p * (exp(tau_p * (SL_curr - SL_slack)) - 1);
else
    F_parallel = -k_p * (exp(tau_p * (SL_slack - SL_curr)) - 1);
end
