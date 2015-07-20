%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: getDeriv                                                  %
%   Date Started: 8/7/2008                                              %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function evaluates the matrix-vector product A*x
%   where A is the matrix of kinetic coefficients and x is the vector
%   of state variables.  This product is the time derivative of state
%   variables, xdot.  Because this system represents the entire model
%   (both iRU and eRU coupled together), the coupling quantities phi
%   and psi must be calculated from x first off.  Then the quantities
%   Ca_i, phi, and psi must be incorporated into the final version of A
%   prior to computing the product A*x.  For efficiency, the iRU model
%   and associated inputs have been condensed to remove Tm state
%   distinctions.  See Program Glossary for variable definitions.
%
%   4/13/09 - This code has been modified to work with the non-XB explicit
%   formulation.  XB binding is computed as a single ODE and has been added
%   to the iRU matrix and vector.
%
%   1/3/11 - Updated with mechanics equations to allow simulation of simple
%   loads, including internal myocyte loads.  Added SL as one of the state
%   vars.
%
%   12/28/11 - Causing the full-length Ca_i to be passed in, rather than
%   determined by a call each time.  linInterp will then be used to extract
%   the proper value.
%***********************************************************************%

function xdot = getDeriv(t,         ... % Current time
                         x,         ... % Vector of state variables
                         AiRU,      ... % Kinetic matrix for iRU model - UNFINISHED
                         Akb,       ... % k_b component of AiRU
                         Agxb,      ... % gxb component of AiRU
                         Ahf,       ... % hf component of AiRU
                         Akon,      ... % k_on component of AeRU
                         Akoff,     ... % k_off component of AeRU
                         alphas,    ... % Fraction of s3 = 1 RU's in each eRU state
                         Ca_i,      ... % Ca transient timecourse
                         t_Ca,      ... % Vector of time points corresponding to Ca_i values
                         L_params,  ... % Parameters for getSL to determine current SL
                         bparams) 
xdot    = zeros(length(x),1); % Initialize vector of time derivatives

[xiRU xeRU] = splitX(x);                    % Split x into iRU and eRU parts    

SL   = xeRU(end - 2);                       % Trim out SL variable
xMpr = xeRU(end - 1);                       % Trim out distortion vars
xMpo = xeRU(end);
xeRU = xeRU(1:end-3);

lenxiRU = length(xiRU);                     % Length of xiRU portion of x                 
B       = sum(xeRU .* (1 - alphas));        % Current occupancy of B state
CM      = sum(xeRU .* alphas);              % Current occupancy of C/M lumped state

xiRU(2) = B - xiRU(1);                      % Update xiRu version of [B1] using xeRu and xiRU info
xiRU(3) = CM - sum(xiRU(4:5));              % Update xiRU version of [C] using xeRU, xiRU, and mass conservation
C       = xiRU(3);
Mpr     = xiRU(4);
Mpo     = xiRU(5);

% Extracting params from bparams
f        = bparams.f;
hb       = bparams.hb;
hf       = bparams.hf;
xMpo0    = bparams.x0;
nu       = bparams.nu;
k_s      = bparams.k_s;
tau_s    = bparams.tau_s;
l_slack  = bparams.l_slack;
k_p      = bparams.k_p;
tau_p    = bparams.tau_p;
SL_slack = bparams.SL_slack;
kxb      = bparams.kxb;


phi     = calcPhi(xiRU); % from state variables
mu      = calcMu(xiRU);


Ca_i  = linInterp(t, t_Ca, Ca_i);       % Obtain current Ca_i
l_tot = Ldep_getLtot(t, L_params);      % Get Current overall length

F_act   = calcActiveForce(x, kxb);      % Calculates the active tension from the filament overlap, distortion variables, and occupancies, scaled by kxb

if ischar(k_s)
    dSL_dt = 0; % This is the constant sarcomere length case
else
    dSL_dt  = (1 / nu) * (calcFSeries(k_s, tau_s, l_slack, l_tot - SL) - calcFParallel(k_p, tau_p, SL_slack, SL) - F_act);   % Compute SL derivative
end

AiRU    = updateAiRU(AiRU,  Akb, Ca_i, Agxb, Ahf, xMpr, xMpo, bparams);  % Update kinetic 
AeRU    = finalizeAeRU(Akon, phi, Akoff, mu);                            % matrices


xdot(1:lenxiRU)         = AiRU * xiRU;       % Calculate derivatives through
xdot(lenxiRU+1:end - 3) = AeRU * xeRU;       % matrix-vector products

xdot(end - 2) = dSL_dt;

xdot(end - 1) = -(f * (C / Mpr) + hb * (Mpo / Mpr)) * xMpr + (dSL_dt / 2);   %xMprdot 
xdot(end)     = -hf * (Mpr / Mpo)*(xMpo - xMpo0) + (dSL_dt / 2);             %xMpodot 


return


