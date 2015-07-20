%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: calcActiveForce                                           %
%   Date Started: 8/26/2008                                             %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function calculates the active force produced based
%   on the assumption that any RU in the 'open' state is producing force.
%   
%   It takes into account filament overlap as a function of sarcomere
%   length according to basic filament geometry.  It also takes into
%   account the distortion of crossbridge populations a-la Ken B.
%   Campbell's work.
%
%   The force is scaled by kxb, which is crossbridge stiffness and the
%   number of XBs all rolled into one tidy scaling factor. It is altered by
%   RLC phosphorylation at the beginning, in getParams.m
%
%   
%***********************************************************************%

function force = calcActiveForce(x, kxb)

[rows, cols] = size(x);

if cols == 1    % If this is a single time point...
    x = x';     % Transpose x to row vector
    num_t_pts = 1;
else
    num_t_pts = cols;
end

[xiRU xeRU]  = splitX(x);

% Setting up vectors
Mpr_vect  = xiRU(:,4);
Mpo_vect  = xiRU(:,5);

SL        = xeRU(:,end-2);
xMpr_vect = xeRU(:,end-1);
xMpo_vect = xeRU(:,end);

% Calculate the overlap
SL_overlap = zeros(num_t_pts, 1);

for j = 1:1:length(SL)
    SL_overlap(j) = getOverlap(SL(j));
end

% Zero out negative distortion values
xMpr_vect = (xMpr_vect > 0) .* xMpr_vect;

% Calculating force
force   = kxb * SL_overlap .* (xMpo_vect .* Mpo_vect + xMpr_vect .* Mpr_vect);

return

