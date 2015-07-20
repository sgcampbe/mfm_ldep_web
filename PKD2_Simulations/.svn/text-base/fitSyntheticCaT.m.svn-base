function [t_fit Ca_fit tau_1 tau_2] = fitSyntheticCaT(CaT_ttp, ...  Should be in ms!
                                                      tau_Ca, ...   Should be in ms!
                                                      Ca_diastolic, ...
                                                      Ca_max, ...
                                                      t_stim)
%
% Function to fit synthetic Ca transients to the mean measured values of
% time to transient peak and the time constant of transient decay.
%
% Uses the synthetic transient of Rice et al 2007.

t = 1:2000;     % Time, in ms.  Assumes the original CaT records are at most 2000 ms in duration

errtol = 0.1;  % Terminate iterations when error less than errtol

% Initial fit values
t_start      = t_stim;
tau_1        = 24;
tau_2        = 159;

init_params  = [tau_1 tau_2]; 

% Assemble fit targets vector
fit_targets = [CaT_ttp, tau_Ca, t_stim];

% Function handle for computing synthetic transient
fhcat = @(tau_1, tau_2) calcCaT(tau_1, tau_2, t, t_start, Ca_diastolic, Ca_max);


% Function handle for error computation
fherr = @(params) evalCaTErr(params, fit_targets, fhcat);


% Call fminsearch
opt = optimset('TolX', errtol, 'TolFun', errtol, 'Display', 'off');
fitted_params = fminsearch(fherr, init_params, opt);


% Recover fitted transient by one more call to calcCaT
tau_1 = fitted_params(1);
tau_2 = fitted_params(2);

[t_fit Ca_fit] = fhcat(tau_1, tau_2);


function err = evalCaTErr(params, fit_targets, fhcat)

% Unpack taus, force to be positive
tau_1 = abs(params(1));
tau_2 = abs(params(2));

stimtype = 'tstim'; % For getCaTProps call

% Unpack targets
CaT_ttp = fit_targets(1);
tau_Ca  = fit_targets(2);
t_stim  = fit_targets(3);

% Compute the synthetic transient
[t Cai] = fhcat(tau_1, tau_2);

% Compute the transient properties
Traw      = t;
Yraw340   = ones(size(t));        % Dummy
Yraw380   = Yraw340 ./ Cai;       % Invert so it looks like a 380 trace
Yrawratio = Yraw340 ./ Yraw380;
props     = getCaTProps(Traw, Yraw340, Yraw380, Yrawratio, 'stimtype', stimtype, 'stim', t_stim);

% Compute the error
err = (tau_Ca - props.tau)^2 + (CaT_ttp - props.t_peak)^2;


fprintf('ttp: %f  tau: %f err: %f \n', props.t_peak, props.tau, err)





function [t Cai] = calcCaT(tau_1, tau_2, t, t_start, Ca_diastolic, Ca_max)

beta = (tau_1 / tau_2)^(-1 / (tau_1 / tau_2 - 1)) - (tau_1 / tau_2)^(-1 / (1 - tau_2 / tau_1));
Cai  = ((Ca_max - Ca_diastolic) / beta) * (exp(-(t - t_start) / tau_1) - exp(-(t - t_start) / tau_2)) + Ca_diastolic;

Cai(t < t_start) = Ca_diastolic;


