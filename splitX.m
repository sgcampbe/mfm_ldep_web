%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: splitX                                                    %
%   Date Started: 8/7/2008                                              %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function assumes a fixed length of xiRU to split the
%   vector of state variables into xiRU and xeRU portions.  Also works
%   on X, the matrix of x as row vectors at time points or Ca points.
%   See Program Glossary for variable definitions.
%
%   7/14/09 - Now this function has been altered to recognize the presence
%   of an additional crossbridge state, meaning that xiRU has 5 entries,
%   not 4.
%***********************************************************************%

function [xiRU xeRU] = splitX(x)

leniRU     = 5;            % UPDATED: now assumes 5 states, B0, B1, C, Mpr and Mpo
[ni nj nk] = size(x);      % Check for timecourse/Ca data

if nj == 1     % x is single column vector
    xiRU    = x(1:leniRU);
    xeRU    = x(leniRU+1:end);   
elseif nk == 1 % Timecourse/Ca data with x in row vectors
    xiRU    = x(:,1:leniRU);
    xeRU    = x(:,leniRU+1:end);
else           % x is 3D with both timecourse and Ca info
    xiRU    = x(:,1:leniRU,:);
    xeRU    = x(:,leniRU+1:end,:);
end


return