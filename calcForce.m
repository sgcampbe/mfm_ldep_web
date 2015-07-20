%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: calcForce                                                 %
%   Date Started: 8/26/2008                                             %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function calculates the total force produced based
%   on the assumption that any RU in the 'open' state is producing force.
%
%   
%***********************************************************************%

function force = calcForce(x, params, time, SL_params)

[rows, cols] = size(x);

[xiRU xeRU]  = splitX(x);


SL         = zeros(rows, 1);
SL_overlap = zeros(rows, 1);

for j = 1:1:length(SL)
    [SL(j) ~] = Ldep_getSL(time(j), SL_params);
    k = getOverlap(SL(j));
    SL_overlap(j) = k;
end


% Setting up vectors
Mpr_vect  = xiRU(:,4);
Mpo_vect  = xiRU(:,5);

SL        = xeRU(:,end-2);
xMpr_vect = xeRU(:,end-1);
xMpo_vect = xeRU(:,end);

% Calculating force
force   = SL_overlap.*(xMpo_vect .* Mpo_vect + xMpr_vect .* Mpr_vect);

% % Plotting SL 
% figure(2)
% plot(time, SL)
% figure(3)
% plot(time, xMpr_vect, 'r')
% figure(4)
% plot(time, xMpo_vect, 'b')
return

