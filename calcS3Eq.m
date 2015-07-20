%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: calcS3Eq                                                  %
%   Date Started: 8/7/2008                                              %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function calculates the total probability
%   associated with the 3 possible s3 (Tropomyosin) states.
%   Handles either column vector xeRU or a matrix of xeRU row vectors
%   at multiple time points (as returned from the ODE solver). Works just
%   as well for multiple Ca points, same format.
%   See Program Glossary for variable definitions.
%
%   4/13/09 - This function has been altered to deal with the simplified
%   models in which the eRU model does not explicitly track XB bound
%   states.  Also, the four state iRU model is assumed.
%
%   7/14/09 - The function has been updated and altered for use in the
%   'linear segment' models, rather than 'ring' models.  I this case, the
%   difference is that myofilament overlap is taken into account. The
%   values in xiRU are only applicable to the population of RUs as a whole
%   for xiRU(1:2) (B0 and B1).  The other states are relative occupancies
%   for the subpopulation of C/M states (they're already normalized to that
%   population's total occupancy).  Thus, the 'real', absolute P{C}, P{M},
%   etc. must come from using this function.
%
%   2/27/10 - Modified to remove lambda effects.
%
%***********************************************************************%

function s3eq = calcS3Eq(x, alphas)

[xiRU xeRU] = splitX(x);

[rows cols] = size(xeRU);           % Check for shape of xeRU

if cols == 1
    mu      = calcMu(xiRU);                   % Fraction of closed/open states that are open
    s3eq(1) = sum(xeRU .* (1 - alphas));        % Probability of blocked RU
    s3eq(2) = mu * (1 - s3eq(1));               % Probability of closed RU
    s3eq(3) = 1 - s3eq(1) - s3eq(2);            % Probability of open RU from conservation of total probability
else
    s3eq    = zeros(rows,3);
    xiRU_raw = xiRU;
    xeRU_raw = xeRU;
    for i = 1:rows
        xiRU      = xiRU_raw(i,:);
        xeRU      = xeRU_raw(i,:)';
        mu        = calcMu(xiRU);                     % Fraction of closed/open states that are open
        s3eq(i,1) = sum(xeRU .* (1 - alphas));        % Probability of blocked RU
        s3gt0     = sum(xeRU .* alphas);              % Probability of NOT blocked RU
        s3eq(i,2) = mu * s3gt0;                       % Probability of closed RU
        s3eq(i,3) = 1 - s3eq(i,1) - s3eq(i,2);        % Probability of open RU from conservation of total probability
    end
end

return

