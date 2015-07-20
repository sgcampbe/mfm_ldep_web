%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: calcMu                                                    %
%   Date Started: 4/14/2009                                             %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function calculates the conditional probability
%   P{s3>0|s3=1}, or mu.  This means the fraction of [C/M] states which are
%   in the closed states, or in other words, the fraction of the lumped
%   state [C/M] which can undergo the transition back to blocked.  This
%   value is needed because XB binding status is no longer tracked by the
%   eRU model.  The rate constant k_off gets multiplied by mu to reflect
%   XB's which are bound and preventing RUs from moving back to blocked.
%
%   3/13/10 - Modified to use xiRU instead of s3eq - ASSUMES 3-STATE XB
%   MODEL IS BEING USED!!!
%***********************************************************************%

function mu = calcMu(xiRU)

mu = xiRU(3) / sum(xiRU(3:5));

return