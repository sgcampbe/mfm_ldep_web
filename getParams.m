%***********************************************************************%
%   Markov Myofilament Model Parameter Setup Function                   %
%   Modified: 7/28/08
%   Author: Stuart Campbell
%   Description: This function sets up all of the free model parameters
%   without doing almost any calculations on them.  ALL FREE PARAMETERS
%   ARE CONTAINED IN THIS FUNCTION.  Exceptions may be passed in through
%   argkey and argvals.
%   Rather than specifying kinetic coefficients directly, this function 
%   may use equilibrium constants and scaling coefficients in some cases.
%   getParams is called only at the highest level (simulation scripts)
%   with one exception:  It is called in solveSingleFullSS to update iRU
%   matrices with Tm kinetic coefficients calculated from the eRU model.
%   In that case, it should be noted that all parameters returned by this
%   function with the exception of kon, koff, f, and g (the Tm parameters,
%   which are actually passed in) are dummy, and are only along for the
%   ride.
%   The output of this function may be used to keep a record of
%   parameter sets which accompany results.
%   9/15/08 - The function has been changed to work with the 4 state iRU
%   scheme.
%
%   4/14/09 - The function has now been further altered to work with the
%   new eRU model which does not explicitly account for XB binding.
%
%   7/14/09 - Altered again for the 'linear' RU arrangement models -
%   includes several new parameters to accommodate a 3-state XB cycling
%   model.  WARNING: SOME XB STRAIN SENSITIVITY PARAMETERS ARE HARD-CODED
%   IN THE updateAiRU FUNCTION!!
%
%   1/03/11 - Added params for simple mechanics model.  Used nonlinear
%   exponential elastic elements a-la Rice 2008.
%
%   3/17/11 - Adding export_params option which will dump all model
%   parameters to 'parameters.txt' in the local directory if set to true.
%
%   5/6/11 - Changed parameterization of RLCP effects according to xb
%   stiffness change.
%***********************************************************************%

function param_out = getParams(varargin)

argkey = {};    % Default, empty argkey
export_params = false;  % Default is false, set this flag to true via 3rd input arg for exporting of parameter values to text file ('parameters.txt' in the calling directory)

% Unpack input arguments
func_argkey = {'argkey', 'argvals', 'export_params'};
for i = 1:nargin
    eval(strcat(func_argkey{i}, ' = varargin{i};'))
end


%-----------------------------------------------%
%---------------Model Parameters----------------%
%-----------------------------------------------%

kb      = 0.09;   %[1/(µM * ms)] - Binding rate of Ca2+ to TnC
ku      = 0.45;   %[1/ms]        - Rate of IpTnI binding to actin and Ca dissociation

kon     = 300.0;   %[1/ms]        - Base rate for transition of Tm from nonpermissive to permissive (k_b+ in the manuscript)
koff    = 0.3270;  %[1/ms]        - Base rate for transition of Tm from permissive to nonpermissive (k_b- in the manuscript)

% Cooperativity parameters
gammaB  = 350;     %[]            - Multiplier enhancing the influence of neighboring RU on activation
q       = 0.5;    %[]            - Weights effects of cooperative coefficients towards forward rates ( q = 1 ) or reverse rates ( q = 0 )

q_Ca    = 1;      %[]            - Number of low affinity Ca binding sites on TnC

Temp    = 273+15; %[K]           - Temperature of experiment

% 3-state Rice XB model parameters
fapp    = 500e-3;     % (1/ms) XB on rate
gapp    = 70e-3;      % (1/ms) XB off rate
hf      = 2000e-3;    % (1/ms) rate between pre-force and force states
hb      = 400e-3;     % (1/ms) rate between pre-force and force states
gxb     = 70e-3;      % (1/ms) ATP consuming transition rate
xbmodsp = 1.0;        % mod to change species 
Q10     = 6.25;       % Temperature dependence of XB kinetics
x0      = 8e-3;       % (µm) strain induced by head rotation
eps_BS  = 2;          % (pN nm^-1) Baseline XB stiffness
kxb     = 10;         % (basically arbitrary force scaling factor)
sigman  = 1;          % Strain dependent parameter - neg strain
sigmap  = 8;          % Strain dependent parameter - pos strain
hfmdc   = 5;          % Strain dependent parameter for forward powerstroke rate

% RLC phoshporylation-related parameters
RLCP     = 0.0;        % Fraction of phosphorylated RLC
p_f      = 0.0;        % Weighting of RLC phos. effect on f
eps_RLCP = eps_BS;     % Stiffness of myosin crossbridge with phosphorylated RLC

% Mechanics/Loading parameters (RFU = relative force units)
nu       = 0.003;       % [RFU * s * um^-1] - Internal viscosity
k_s      = 'isometric'; % [RFU * um^-1]     - Series elasticity - set equal to the string 'isometric' for isometric contraction; set to zero for unloaded shortening; in between to mimic end compliance of prep
tau_s    = 5;           % [um^-1]           - Exponential length dependence constant for nonlinear series stiffness (default value here is arbitrary)
l_slack  = 2;           % [um]              - The absolute value of this parameter shouldn't be very critical - it just determines where slack is for the series elastance
k_p      = 0.002;       % [RFU * um^-1]     - Parallel (internal) elasticity, basically titin stiffness (from Rice 2008)
tau_p    = 10;          % [um^-1]           - Exponential length dependence constant for nonlinear parallel stiffness (from Rice 2008)
SL_slack = 1.9;         % [um]              - Slack sarcomere length, as chosen by Rice et al. 2008
kxb      = 100;         % []                - Lumped active tension scaling factor, representing crossbridges and the number of cross bridges.  Should generally be adjusted until the max active response is ~1


%--------------------------%
% Unpack Passed Parameters %
%--------------------------%
% This replaces any default parameters defined above with ones passed into
% the function.
if nargin
    for i = 1:length(argkey)
        eval(strcat(argkey{i}, '= argvals{i};'))
    end
end 

%--------------------------------%
% Calculate rate pairs using K's %
%--------------------------------%
% If equilibrium constants are passed in, update 'off' rates
if isInList('Kb', argkey) 
    ku     = kb     / Kb;                                   %[1/ms] - Unbinding rate of Ca2+ from TnC
end
if isInList('Kon', argkey)
    koff   = kon    / Kon;                                  %[1/ms] - Base rate for transition of Tm from permissive to nonpermissive
end

%-----------------------------------%
% Other Calculations for XB cycling %
%-----------------------------------%

Q = xbmodsp * Q10^((Temp-310)/10);  % Lumped all-purpose scaling factor for temperature and species...

fappT   = fapp * Q;
gappT   = gapp * Q;
hfT     = hf   * Q; % Next 3 get scaled by distortion-dependent processes later
hbT     = hb   * Q;
gxbT    = gxb  * Q;

%------------------%
% Add RLCP Effects %
%------------------%

fappT_raw   = fappT;    % Save Copy for export 
hfT_raw     = hfT;
hbT_raw     = hbT;

fappT    = fappT_raw * (1 + p_f  * RLCP);
eps_avg  = eps_BS * (1 - RLCP) + eps_RLCP * RLCP;    % Weighted average of stiffness
kinetics = computePowerstrokeMechanoEnergetics(eps_avg, hfT, hbT, eps_BS, x0 * 1e3, Temp);  % Note conversion of x0 to nm in ARG LIST!!
hfT      = kinetics.h_f_RLCP;
hbT      = kinetics.h_b_RLCP;


%-------------------------------%
% Update kxb to reflect eps_avg %
%-------------------------------%

kxb = (kxb / eps_BS) * eps_avg; % kxb is implicitly composed of the number of myosin XB's per c.s. area TIMES the crossbridge stiffness


%----------------------------%
% Package params for passing %
%----------------------------%
% Package params for passing
% NOTE THAT THIS ORDER OF ASSIGNMENT MUST AGREE WITH THAT USED IN THE
% FUNCTION getiRUParams FOR THE FIRST 9 VARS!!!

% Ca binding
param_out.kb    = kb;
param_out.ku    = ku;

% Tropomyosin shifting
param_out.kon   = kon;
param_out.koff  = koff;

% XB cycling
param_out.f     = fappT;
param_out.g     = gappT;
param_out.hf    = hfT;
param_out.hb    = hbT;
param_out.gxb   = gxbT;

% This packing order is not critical
param_out.q_Ca    = q_Ca;
param_out.gammaB  = gammaB;
param_out.q       = q;
param_out.x0      = x0;
param_out.nu      = nu;
param_out.k_s     = k_s;
param_out.tau_s   = tau_s;
param_out.l_slack = l_slack;
param_out.k_p     = k_p;
param_out.tau_p   = tau_p;
param_out.SL_slack= SL_slack;
param_out.kxb     = kxb;
param_out.eps_avg = eps_avg;
param_out.sigman  = sigman;
param_out.sigmap  = sigmap;
param_out.hfmdc   = hfmdc;

% Export parameters, if selected
if export_params
    % Dump some intermediate variables too 
	param_out.f_raw    = fappT_raw;
    param_out.hf_raw   = hfT_raw;
    param_out.hb_raw   = hbT_raw;
    param_out.p_f      = p_f;
    param_out.eps_BS   = eps_BS;
    param_out.eps_RLCP = eps_RLCP;
    param_out.T        = Temp;
    
    % VERY most raw XB values
    param_out.f_base   = fapp;
    param_out.g_base   = gapp;
    param_out.hf_base  = hf;
    param_out.hb_base  = hb;
    param_out.gxb_base = gxb;
    param_out.xbmodsp  = xbmodsp;
    param_out.RLCP     = RLCP;
    
    %Precompute some things
    param_out.koffBB_raw = koff * (gammaB^-2)          ^ -(1 - q);
    param_out.koffBC_raw = koff * (gammaB^-1)          ^ -(1 - q);
    param_out.koffCC_raw = koff * (1)                  ^ -(1 - q);
    param_out.konBB_raw  = kon  * (gammaB^-2)          ^ q;
    param_out.konBC_raw  = kon  * (gammaB^-1)          ^ q;
    param_out.konCC_raw  = kon  * (1)                  ^ q;
    
       
    params_to_write = rmfield(param_out, 'k_s');  % It's not needed in Continuity, and can mess up output...
    
    % Write to text file
    pnames = fieldnames(params_to_write);
    fptr = fopen('parameters.txt', 'w');
    for i = 1:length(pnames)
        fprintf(fptr, '%s\t%.13E\n', pnames{i}, params_to_write.(pnames{i}));
    end
    fclose(fptr);
end

return