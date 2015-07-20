%************************************************************************%
%   Length-dependent Myofilament Model
%   File:   Ldep_getLtot.m                                               %
%   Date Started: 1/3/2011                                               %
%   Author: Stuart Campbell
%   Description: This function gives the total length of the simulated
%   muscle, based on the time, t, and the parameters in L_params.
%*********************************************************************************%


function L = Ldep_getLtot(t, L_params)

L_code = L_params{1};     
Lmax   = L_params{2};     %Max L that is reached during stretch
Lmin   = L_params{3};     %Length of iso segment before PS
iso_t1 = L_params{4};     %Length of iso segment before PS
strDur = L_params{5};     %Duration of stretch in ms
        
switch L_code
    case 1      % Isometric about Lmin
        L = Lmin;
        
    case 2      % "after" prestretch MuscleLength data taken from 129 mouse (100329_m02 - PS_5pct_after_PS) and altered to reach L = 2.45 
        
        if (t < iso_t1)               %Not sure if MATLAB syntax allows this in an if statement
            L     = Lmin;
        elseif(t > strDur + iso_t1)
            L     = Lmin;
        else
            L = (Lmax - Lmin) * 0.5 * (1 - cos(2 * (pi / strDur) * (t - iso_t1))) + Lmin;
        end
                     
end

return