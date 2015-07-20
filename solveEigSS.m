%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: solveEigSS                                                %
%   Date Started: 8/26/2008                                             %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function essentially saves one line of code each
%   time a steady-state solution from solving an eigenvalue problem is
%   desired.  Because the corresponding eigenvector is not scaled, it
%   must be nomalized such that its sum is equal to one prior to being
%   used in any of the sumulation contexts.
%***********************************************************************%

function x0 = solveEigSS(A)

opts.disp = 0;                      % Mute eigs output

[xraw junk] = eigs(A,1,1e-10,opts); % Solve eigenvalue problem
    
x0 = xraw ./ sum(xraw);             % Normalize so sum to 1

return