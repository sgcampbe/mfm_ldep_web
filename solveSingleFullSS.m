%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: solveSingleFullSS                                         %
%   Date Started: 8/25/2008                                             %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function uses a series of eigenvalue problem
%   solutions to calculate the steady-state values of the vector of
%   state variables, xInit, for the full iRU/eRU model.  It does this
%   irrespective of any pruning.  In STEP 1, iRU model is solved to 
%   determine phi at the given level of Ca.  In STEP 2, the eRU model is 
%   solved for that phi.  In STEP 3, kon, koff, f, and g are
%   updated via equilibrium values determined from solution of the eRU
%   model at constant phi, and the iRU model is solved once again to
%   obtain the final values.  While phi will be the same after the second
%   solve, the values of the state variables will have been updated to
%   their correct steady-state values.
%   This function may be conveniently used to generate initial
%   conditions for dynamic simulations.
%
%   4/14/09 - Small changes made to suit simplified eRU model as well as
%   4-state iRU model.
%***********************************************************************%

function xInit = solveSingleFullSS(params, ...
                                   Ca_i , ...     Ca concentration
                                   SL)         % SL of steady-state solution
                                   
%-------------------------------------%
% Extract quantities from params struct
names = fieldnames(params);
for i = 1:length(names)
    eval(strcat(names{i},'=','params.',names{i},';'));
end
% Now all param names are defined
%-------------------------------------%

x0  = bparams.x0;
xMpr = 0.0;  % At steady state, pre-isomerization XBs have zero distortion
xMpo = x0;   % At steady state, post-isomerization XBs have the default, powerstroke distortion

%----------------------%
% Solve model - STEP 1 %
%----------------------%

xiRU   = solveSSiRU(AiRU, ...                       % Solve steady-state
                    namesp, ...                     % Note that this produces xiRU based on
                    Ca_i, ...                       % neighbor-independent kon/koff, so abs vals are not correct
                    xMpr, ...
                    xMpo, ...
                    x0);
                
phi    = calcPhiSS(xiRU);                           % Calculate conditional probablility phi
mu     = calcMu(xiRU);                              % Calculate conditional probablility mu 
xiRU(3:5) = xiRU(3:5) / sum(xiRU(3:5));             % 'Fix' xiRU to be the coupled version (C, MPr and MPo states are normalized)


%----------------------%
% Solve model - STEP 2 %
%----------------------%

xeRU    = solveSSeRU(Akon, phi, Akoff, mu); % Solve steady-state

s3eq    = calcS3Eq([xiRU; xeRU], alphas);

%----------------------%
% Solve model - STEP 3 %
%----------------------%

xiRU(2) = phi * s3eq(1);        % Prob of B1 state
xiRU(1) = s3eq(1) - xiRU(2);    % Prob of B0 state
xiRU(3:5) = xiRU(3:5) * (1 - s3eq(1));             % 'Fix' xiRU to be the coupled version (C, MPr and MPo states are normalized)

%---------------------%
% Return final result %
%---------------------%

xInit = [xiRU; xeRU];                              % Join iRU and eRU results to form single vector of st var values

return
