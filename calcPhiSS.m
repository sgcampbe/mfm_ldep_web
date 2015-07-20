%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: calcPhiSS                                                 %
%   Date Started: 8/25/2008                                             %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function uses iRU state variable values to
%   calculate the quantity phi.  The function calcPhi uses a combination
%   of iRU and eRU information to do the same.  iRUmap, a list of
%   occupied or permitted iRU states, is used to account for the fact that
%   the iRU model may be pruned to some degree to simulate some type of
%   experimental condition.
%   This function can handle either a single column vector of state vars
%   or alternatively a matrix where rows are state var vecs at different
%   concentrations of Ca.
%
%   4/13/09 - Changed to work exclusively with 4 state iRU model.
%***********************************************************************%

function phi = calcPhiSS(x)

[rows cols] = size(x);      % Get dimensions of x

if cols > 1                 % If x contains multiple Ca points
    x = x';                 % Transpose x 
    [rows cols] = size(x);  % Update rows/cols
    phi = zeros(cols,1);    % Dimension phi output vector
end

for i = 1:cols
    phi(i) = x(2,i) ...     % Sum of blocked RUs w/ Ca bound to
           / sum(x(1:2,i)); % Divided by sum of all blocked RUs
end

