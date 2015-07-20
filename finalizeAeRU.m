%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: finalizeAeRU                                              %
%   Date Started: 8/21/2008                                             %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function accepts the kinetic matrix Akon and the  
%   partially-assembled matrix AeRU (containing the sum of the other three matrices)
%   and the quantity phi and returns the finished matrix.
%   See Program Glossary for variable definitions.
%
%   4/14/09 - Modified for the simpler eRU model which does not represent
%   XB binding explicitly.  
%***********************************************************************%

function AeRU = finalizeAeRU(Akon, phi, Akoff, mu)

AeRU = phi * Akon + mu * Akoff;

return
