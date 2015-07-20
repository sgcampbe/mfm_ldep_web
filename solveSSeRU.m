%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: solveSSeRU                                                %
%   Date Started: 8/26/2008                                             %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function solves the eRU model for given values of phi
%   and mu (updated for non-explicit XB binding eRU model, 4/14/09).
%
%***********************************************************************%

function xeRU = solveSSeRU(Akon, phi, Akoff, mu)

AeRU      = finalizeAeRU(Akon, phi, Akoff, mu); % Update eRU matrix to current phi
xeRU      = solveEigSS(AeRU);                   % Solve steady-state

return