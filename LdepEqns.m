% Differential Equations for Ldep model

function ydot = LdepEqns(t,         ... Current time
                         y,         ... Vector of state variables at the current time
                         params,    ... Cell array of model parameter values
                         Ca_params, ... Cell array of Calcium transient parameters
                         SL_params)   % Cell array of SL pattern parameters

ydot = zeros(length(y), 1);                % Set ydot up as a column vector (required by ode23t)

Ca = Ldep_getCa(t, Ca_params);            % Get Current Ca concentration
[SL dSL_dt] = Ldep_getSL(t, SL_params);   % Get Current SL



%-------------Set up model parameters--------------------%
gammaB = params.gammaB;
k_negB = params.k_negb;
k_posB = params.k_posb;

k_negCa = params.k_negCa;
k_posCa = params.k_posCa;

f = params.f;
g = params.g;

gxb0 = params.gxb0;

hf = params.hf;
hb = params.hb;

xpo0 = params.xpo0;

sigmap = params.sigmap;
sigmam = params.sigmam;

%Grouping all state variables into single array y
eru1 = y(1);
eru2 = y(2);
eru3 = y(3);
eru4 = y(4);
B0   = y(5);
Mpr  = y(6);
Mpo  = y(7);
xpr  = y(8);
xpo  = y(9);

% Zero out xpr and xpo if negative (for purposes of calculations here only)
% They have to be zeroed out in the driving script after the simulation
% also for force to be correct.
if xpr < 0
    xpr = 0;
end

if xpo < 0
    xp0 = 0;
end


% Update distortion dependence
sigma = (xpo < xpo0) * sigmam + (xpo <= xpo0) * sigmap; % Update sigma depending on mean distortion (positive or negative)
gxb = gxb0*exp(sigma*((xpr - xpo0)^2)); %eqn 13 from Razumova paper, makes the rate of cross bridge detachment strain dependent 

%----------Algebraic equations (determine values for 'missing' states)---%
B   = [eru1 eru2 eru3 eru4] * [1 2/3 1/3 0]';      % Current occupancy of B state - note the fact that I'm using a dot product here
CM  = [eru1 eru2 eru3 eru4] * [0 1/3 2/3 1]';      % Current occupancy of C/M lumped state
B1  = B - B0;                                      % Conservation of mass (conservation of blocked states)
C   = CM - (Mpr + Mpo);                            % Conservation of mass (for lumped CM state)
psi = C / CM;                                      % Fraction of CM states that are crossbridge free (able to return to blocked state)
phi = B1 / B;                                      % Fraction of blocked states that have Ca bound and can therefore advance to Closed/Open state


%-------------------Calculations of derivatives--------------------------%
deru1_dt = psi*gammaB*k_negB*eru2 - 3*phi*(1/gammaB)*k_posB*eru1;
deru2_dt = 3*phi*(1/gammaB)*k_posB*eru1 + 2*psi*sqrt(gammaB)*k_negB*eru3 - 2*phi*sqrt(gammaB)*k_posB*eru2 - psi*gammaB*k_negB*eru2;
deru3_dt = 2*phi*sqrt(gammaB)*k_posB*eru2 + 3*psi*k_negB*eru4 - 2*psi*sqrt(gammaB)*k_negB*eru3 - phi*k_posB*eru3;
deru4_dt = phi*k_posB*eru3 - 3*psi*k_negB*eru4;

dB0_dt = k_negCa*B1 - k_posCa*Ca*B0;    % Stuart added Ca concentration to the second term
dMpr_dt = f*C + hb*Mpo - hf*Mpr - g*Mpr;
dMpo_dt = hf*Mpr - hb*Mpo - gxb*Mpo;

dxpo_dt = -hf*(Mpr/Mpo)*(xpo - xpo0) + (dSL_dt / 2);         %eqn 11 from Razumova paper - Note that dSL_dt is cut in half to reflect half-sarcomere length
dxpr_dt = -(f*(C/Mpr) + hb*(Mpo/Mpr))*xpr + (dSL_dt / 2);    %eqn 12 from Razumova paper


%Packing calculated derivatives into vector ydot
ydot(1) = deru1_dt;
ydot(2) = deru2_dt;
ydot(3) = deru3_dt;
ydot(4) = deru4_dt;
ydot(5) = dB0_dt;
ydot(6) = dMpr_dt;
ydot(7) = dMpo_dt;
ydot(8) = dxpr_dt;
ydot(9) = dxpo_dt;












