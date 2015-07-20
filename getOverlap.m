%***********************************************************************%
% Function: getLambdas.m
% Author: Stuart Campbell
% Date Started: 7/14/09
% Description: Returns the fraction of single overlap for an average thin
% filament.  The force is assumed to be proportional to the probability of
% XB binding scaled by the overlap.
%
%***********************************************************************%

function ol = getOverlap(SL)

%------------%
% Parameters %
%------------%

N = 26; % Number of RUs per thin filament

% Mainly from Rice, 2008 model
ltn    = 1.2;  % [um] - length of thin filament from z to m line
ltk    = 1.65; % [um] - length of thick filament, from tip to tip
lbare  = 0.1;  % [um] - width of bare zone on thick filament

sigma   = 50;    % Parameter determining steepness of sigmoid functions at boundaries

%--------------%
% Calculations %
%--------------%

hSL = 0.5 * SL;     % Pre-compute half SL

% Determine half-peak value of z-line end single overlap (l50ze)
l50ze = max(0, hSL - (ltk / 2));

% Determine half-peak value of m-line end single overlap (l50me)
l50me_raw = [ltn, ...            % For when hSL > ltn + 0.5 * lbare
             hSL - lbare/2, ...  % For when ltn-lbare/2 < hSL < ltn + 0.5 * lbare
             SL - ltn];          % For when hSL < ltn-lbare/2
l50me = min(l50me_raw);          

% Get lambdas
n_ltns  = (((1:N) - 0.5) ./ N) * ltn;  % Determine thin filament position for each RU
lambdas = makeSig(n_ltns, l50ze, sigma).*makeSig(n_ltns, l50me, sigma, -1, 1);

ol = sum(lambdas) / N;

return