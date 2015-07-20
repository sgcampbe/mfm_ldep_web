function kinetics = computePowerstrokeMechanoEnergetics(eps_RLCP, h_f_BS, h_b_BS, eps_BS, x0_BS, T)

%***********************************************************************%
%   Computing Statistical Mechanics of RLC Phosphorylation Mechanism    %
%   Script: script_computePowerstrokeMechanoEnergetics.m                %
%   Date Started: 5/5/2011                                              %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function performs computations to determine how
%   crossbridge cycling rates would be affected by lever arm stiffening due
%   to regulatory light chain phosphorylation.  It can be called by the
%   parameter function (getParams.m) in order to update h_f and h_b for
%   given RLC-phosphorylation induced changes to the XB stiffness epsilon.
%***********************************************************************%

% Description of input parameters
% h_f_BS  [s^-1]     - Forward pwr stroke rate
% h_b_BS  [s^-1]     - Reverse pwr stroke rate
% eps_BS  [pN*nm^-1] - myosin XB stiffness
% x0_BS   [nm]       - Pwr stroke distance
% T       [K]        - Temperature

% Constants
R  = 8.314472;      % [J K^1 mol^1] - Gas Constant
kB = 1.3806504e-23; % [J K^-1]        - Boltzmann Constant
h  = 6.6260690e-34; % [J s]           - Planck Constant

% Convert kinetic constants to seconds
h_f_BS = h_f_BS * 1000;
h_b_BS = h_b_BS * 1000;


% Calculate pre-post pwr stroke energy difference, baseline stiffness
% Via Gibbs relation
K_prpo_BS  = h_f_BS / h_b_BS;           % Equilibrium constant
dG_prpo_BS = -log(K_prpo_BS) * R * T;   % Change in state energy via Gibbs

% Calculate zero-stiffness (zs) energy difference between pre-post states,
% by assuming linear XB stiffness
dG_spr_BS = calcdG_spr(eps_BS, x0_BS);      % [J mol^-1] - Compute work per mol of XB's

dG_prpo_ZS = dG_prpo_BS - dG_spr_BS;        % [J mol^-1] - zero stiffness energy difference

% Calculate forward activation energy based on forward rate (Eyring Eqn)
dG_f_act_BS = -log((h_f_BS * h) / (kB * T)) * R * T;

% Calculate zero-stiffness forward activation energy assuming that peak
% energy barrier under zero stiffness coincides with half of the
% powerstroke distance
dG_f_act_ZS = dG_f_act_BS - calcdG_spr(eps_BS, 0.5 * x0_BS);    % NOTE: FIXED 1/2 mistake found by reviewer!

% % Calculate zero-stiffness reverse activation energy 
dG_b_act_ZS = dG_f_act_ZS - dG_prpo_ZS;


% Now, it is possible to calculate the kinetic rates as a function of the
% crossbridge stiffness, eps_XB.  The first step is to update the
% activation energies according to the XB strain energy, then the second 
% is to apply the Eyring equation

dG_f_act = dG_f_act_ZS + calcdG_spr(eps_RLCP, 0.5 * x0_BS);                                 % NOTE: FIXED 1/2 mistake found by reviewer!
dG_b_act = dG_b_act_ZS - calcdG_spr(eps_RLCP, x0_BS) + calcdG_spr(eps_RLCP, 0.5 * x0_BS);   % NOTE: FIXED a different bug having to do with calculation of the reverse rate!


h_f = ((kB * T) / h) * exp(-dG_f_act / (R * T));
h_b = ((kB * T) / h) * exp(-dG_b_act / (R * T));

% Calculate pre-post pwr stroke energy difference, PHOSPHORYLATED stiffness
% Via Gibbs relation - THIS CALC IS TO BE USED IN THE MANUSCRIPT!
K_prpo_PHOS  = h_f / h_b;                   % Equilibrium constant
dG_prpo_PHOS = -log(K_prpo_PHOS) * R * T;   % Change in state energy via Gibbs
fprintf('Forward powerstroke activation energy                    : %f\n', dG_f_act)
fprintf('Relative energy change, pre to post powerstroke, non-phos: %f\n', dG_prpo_BS)
fprintf('Relative energy change, pre to post powerstroke, PHOS    : %f\n', dG_prpo_PHOS)
fprintf('Phosphorylation increases post-powerstroke energy by     : %f\n', dG_prpo_PHOS - dG_prpo_BS)

% Store kinetic results
kinetics.h_f_RLCP = h_f / 1000; % NOTE: Convert back to 1/ms!!!
kinetics.h_b_RLCP = h_b / 1000; % NOTE: Convert back to 1/ms!!!

return

% Helper function to calculate potential energy per mol of post-powerstroke
% XBs.  NOTE: stiffness of XB (eps_XB) should be in pN nm^-1, head rotation
% (x0) should be in nm.
function dG_spr = calcdG_spr(eps_XB, x0)

NA = 6.02214179e23; % [mol^-1] - Avogadro Constant

ps_work   = (1 / 2) * eps_XB * (x0 ^ 2);    % [pN * nm] - Work by one powerstroke 
ps_work_J = ps_work * (1e-12) * (1e-9);     % [J]  - Convert to joules
dG_spr    = ps_work_J * NA;                 % [J mol^-1] - Compute work per mol of XB's