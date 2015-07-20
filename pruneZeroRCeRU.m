%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: pruneZeroRCeRU                                            %
%   Date Started: 8/21/2008                                             %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function accepts the four separate kinetic matrices 
%   for the AeRU model and returns new versions in which any rows or 
%   columns containing entirely zeros are eliminated, and the size of 
%   the matrix reduced.
%   See Program Glossary for variable definitions.
%
%   4/14/09 - Modified for the simpler eRU model which does not represent
%   XB binding explicitly.  Now only Akon and Akoff are inputs.
%***********************************************************************%

function [Akon Akoff] = pruneZeroRCeRU(Akon, Akoff)

thresh = 1e-8;      % Artificial threshold for roundoff error

AeRU = Akoff;
AeRU = finalizeAeRU(Akon, 1.0, Akoff, 1.0);

rows = sum(abs(AeRU), 2) > thresh;
cols = sum(abs(AeRU), 1) > thresh;

Akon  = Akon(rows, cols);
Akoff = Akoff(rows, cols);

return
