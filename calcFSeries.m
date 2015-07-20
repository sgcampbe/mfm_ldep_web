%************************************************************************%
% Myofilament Model - Length Dependent
% Function: calcFSeries.m
% Author: Stuart Campbell
% Start date: 12/29/10
% Description: Calculate the force due to series elastance after the
% formulation of Rice (in terms of the exponential formulation)
%************************************************************************%

function F_series = calcFSeries(k_s, tau_s, L_slack, L_curr)

if ischar(k_s)
    error('calcFSeries was called with k_s = %s, the calling function must catch this!', k_s)
else
    if L_curr >= L_slack
        F_series = k_s * (exp(tau_s * (L_curr - L_slack)) - 1);
    else
        F_series = k_s * (exp(tau_s * (L_curr - L_slack)) - 1);
        %F_series = 0;
    end
end