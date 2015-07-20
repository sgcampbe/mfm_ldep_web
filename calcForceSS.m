%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: calcForceSS.m                                             %
%   Date Started: 8/26/2008                                             %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function calculates the total force based on
%   occupancy of the open state, along with the distribution of XBs in that
%   state (pre-rotation vs. post-rotation), and the average distortion of
%   either population (xPr and xPo).  This is all done according to the
%   formulation of Rice et al, 2008, with some modifications.
%
%   This version is for steady-state, when xPo may be assumed to be equal
%   to x0 and the SL is constant.  
%***********************************************************************%

function force = calcForceSS(x,...      % The Markov vector (no extras state vars) 
                             params,... % A full version of the parameters 
                             SL)        % The steady-state SL

kxb     = params.bparams.kxb;
alphas  = params.alphas;
numRUs  = params.numRUs;

x0      = params.bparams.x0;
xPo     = x0;   % Steady state dictates these two lines
xPr     = 0.0;

[rows, cols] = size(x);
[xiRU xeRU]  = splitX(x);

if cols == 1
    MPr     = xiRU(4);
    MPo     = xiRU(5);
    % force   = kxb * (xPo * MPo + xPr * MPr);
    force   = MPo * x0;
else    % must be data at multiple pCa's
    force = zeros(rows,1);
    for i = 1:rows
        MPr      = xiRU(i,4);
        MPo      = xiRU(i,5);
        % force(i) = kxb * (xPo * MPo + xPr * MPr);
        force(i) = MPo * x0;
    end        
end
return

