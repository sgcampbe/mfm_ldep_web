%***********************************************************************%
%   Markov Myofilament Model: Calculate Eigenvalue prob input params    %
%   File:   geteRUParams                                                %
%   Date Started: 8/25/2008                                             %
%   Author: Stuart Campbell                                             %
%   Description: This function uses passed in basic params in conjunction 
%   with getEqnDataX.m (where X is the number of RUs in a ring of the eRU 
%   model) to compute eRU kinetic coefficient matrices. It prunes states
%   according to mdex, with mdex = -1 being the XB-free case. The
%   structures alphas and betas are pruned and placed, ready to go, in the
%   output params structure.
%
%   9/18/08 - Changed output so that bparams struct is added to param_out.
%   This allows one of the alternate ktr schemes to re-calculate a portion
%   of the AeRU matrix.  Also added mdex.
%
%   3/16/09 - Removed mdex-dependent transType lookup.  This doesn't work
%   conceptually, as the model doesn't track individual RUs.
%***********************************************************************%

function param_out = geteRUParams(varargin)

mdex     = 1;   % Default mdex - corresponds to full Markov model

% Unpack input arguments
func_argkey = {'numRUs', 'bparams', 'mdex'};
for i = 1:nargin
    eval(strcat(func_argkey{i}, ' = varargin{i};'))
end

param = getEqnData(numRUs);  % Get dynamically generated structures for ensemble RU model

kon    = bparams.kon;
koff   = bparams.koff;
gammaB = bparams.gammaB;
q      = bparams.q;

nS     = length(param.alphas);                   % Number of ensemble RU model states
Akon   = sparse([],[],[],nS,nS,0);               % Initialize kinetic coefficient matrices
Akoff  = sparse([],[],[],nS,nS,0);               % Initialize kinetic coefficient matrices

% Apply constraints formulation to calculate final kinetic coefficients

gammaM = 1; % gammaM is just dummy, it's not in this formulation

koffBB = koff * (gammaB^-2)          ^ -(1 - q);
koffBC = koff * (gammaB^-1)          ^ -(1 - q);
koffBM = koff * (gammaB^-1 * gammaM) ^ -(1 - q);
koffCC = koff * (1)                  ^ -(1 - q);
koffCM = koff * (gammaM)             ^ -(1 - q);
koffMM = koff * (gammaM^2)           ^ -(1 - q);

konBB  = kon * (gammaB^-2)          ^ q;
konBC  = kon * (gammaB^-1)          ^ q;
konBM  = kon * (gammaB^-1 * gammaM) ^ q;
konCC  = kon * (1)                  ^ q;
konCM  = kon * (gammaM)             ^ q;
konMM  = kon * (gammaM^2)           ^ q;


% Create static look-up table for kinetic parameters based on transType
k_coeff_on   = {'konBB', 'konBC', 'konCC', 'konBM', 'konCM', 'empty', 'konMM'};    
k_coeff_off  = {'koffBB', 'koffBC', 'koffCC', 'koffBM', 'koffCM', 'empty', 'koffMM'};  

% Determine prunelist using mdex
if mdex == -1                               % If this is an XB-free simulation
    prunelist = find(param.betas > 0);      % Eliminates all states with XBs
    mdex = 1;                               % Set mdex to 1 for kinetics determination
    wasneg1 = 1;                            % Sets alternate flag (for alpha/beta pruning)
else
    prunelist = [];
    wasneg1 = 0;
end

%------------------------------%
% Generate kinetic rate matrix %
%------------------------------%

for j = 1:length(param.transTypes)
    ipre  = param.iPre(j);
    ipost = param.iPost(j);
    if (sum(ipre == prunelist) > 0) | (sum(ipost == prunelist) > 0) % State pruning
        continue
    else
        Akon(ipre,  ipre)   = Akon(ipre,  ipre) - param.m(j) * eval(k_coeff_on{param.transTypes(j)});
        Akon(ipost, ipre)   = Akon(ipost, ipre) + param.m(j) * eval(k_coeff_on{param.transTypes(j)});
        Akoff(ipost, ipost) = Akoff(ipost, ipost) - param.n(j) * eval(k_coeff_off{param.transTypes(j)});
        Akoff(ipre,  ipost) = Akoff(ipre,  ipost) + param.n(j) * eval(k_coeff_off{param.transTypes(j)});
    end
end

[Akon Akoff]  = pruneZeroRCeRU(Akon, Akoff);     % Prune out rows/cols corresponding to pruned states

%---------------------------%       (This case must by solved for in a
% Check for 'all s3>0' mdex %       special way if the simulation is 
%---------------------------%       dynamic)

if param.alphas(mdex) > 0.99   % If all RU the mask are closed or open
    isallon = true;
else
    isallon = false;
end

%------------------------%
% Prune alphas and betas %
%------------------------%
% (Pruning is skipped if mdex == 1)

map = logical(ones(nS,1)); %#ok<LOGL>
if (mdex ~= 1) || wasneg1
    map(prunelist) = 0;
end

alphas = param.alphas(map);

%-------------------------------%
% Package structures for return %
%-------------------------------%

param_out.eRUprunelist = prunelist;
param_out.Akon         = Akon;
param_out.Akoff        = Akoff;
param_out.alphas       = alphas;
param_out.numRUs       = numRUs;
param_out.isallon      = isallon;
param_out.bparams      = bparams;
param_out.mdex         = mdex;
param_out.mask_alpha   = param.alphas(mdex); % CAREFUL! Use pre-pruning alphas!!
return