%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: runSSFull                                                 %
%   Date Started: 8/25/2008                                             %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function simulates a steady-state pCa
%   relation using the coupled iRU/eRU model. Raw state variable values
%   are returned from this function for a variety of post-processing
%   functions to operate on (it's general).
%
%   7/15/09 - Modified to conform to linear chain of RUs, 5-state iRU
%   model, and non-isometric formulation.
%***********************************************************************%

function data = runSSFull(params, Ca_min, Ca_max, numPts, SL)

%-------------------------%
% Setup output structures %
%-------------------------%

niRUst = 5;                                 % Number of iRU model states
neRUst = length(params.alphas);             % Number of eRU model states

X      = zeros(numPts, (niRUst + neRUst));  % Initialize matrix of solved st var vectors

%-----------------------%
% Setup Ca Range vector %
%-----------------------%

Ca_range = makeLogCaRange(Ca_min, Ca_max, numPts);

%---------------------%
% Loop over Ca values %
%---------------------%

for i = 1:numPts
    X(i,:) = (solveSingleFullSS(params, Ca_range(i), SL))';
end

%-----------------%
% Package Outputs %
%-----------------%

data.X = X;
data.Ca_range = Ca_range;

return



