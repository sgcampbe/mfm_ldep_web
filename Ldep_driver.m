%**********************************************************************%
%   Length dependent myofilament Model - driver                        %
%   File:   Ldep_driver.m                                              %
%   Author: Jared Tangney                                              %
%   Description: This file will set up the ODEs with intial conditions %      
%   and then calls an ODE solver to solve the ODEs.                    % 
%**********************************************************************%

clear all, close all, clc

%-----------------------------------------------%
%-------------Simulation Parameters-------------%
%-----------------------------------------------%

t_start = 0;    % Start time of simulation (ms)
t_end   = 500;  % End time of simulation (ms)
max_step = 1;   % Maximum integrator step (ms)

%------------------%
% Model parameters %
%------------------%

params.gammaB = 300;
params.k_negb = 0.327;
params.k_posb = 350;

params.k_negCa = 0.45;
params.k_posCa = 0.09;


params.f = 2.22e-2;
params.g = 7.8e-3;

params.gxb0 = 0.22;       %value of gxb at isometric condition 

params.hf = 0.22;
params.hb = 4.44e-2;

params.xpo0 = 0.007;        % [microns] isometric value of average distortion among bridges in post power stroke state 

params.sigmap = 8;          %paramater that grades distortional dependence (sigma = 0 means no distortional dependence)
params.sigmam = 1;          %same parameter, but for shortening conditions (when xpo

Ca_params = {2, 0.1, 1.0};  % Ca_code, Ca_min, Ca_max
SL_params = {2 2.25 2.0 150 150}; % SL_code SLmax SLmin iso_t1 strDur

%--------------------%
% INITIAL CONDITIONS %
%--------------------%

eru1 = 0.8295;
eru2 = 0.0007;
eru3 = 0.0179;
eru4 = 0.1519;

B0  = 0.8196;            

Mpr  = 0.0157;
Mpo  = 0.0131;
xpr  = 0;
xpo  = params.xpo0;

y0 = [eru1 eru2 eru3 eru4 B0 Mpr Mpo xpr xpo];  % package IC's into y0 vector

t_span  = t_start:t_end;                        % Vector containing desired timepoints of solution
options = odeset('MaxStep',max_step,'Stats','on');  % Set integrator options

[T,y]   = ode23t(@(t,y) LdepEqns(t,         ... % Current time
                                 y,         ... % Vector of state variables
                                 params,    ... % The model parameters
                                 Ca_params, ... % Parameters for getCa to determine current Ca_i
                                 SL_params),... % Parameters for getSL to determine current SL
                 t_span, y0, options);

%-----------------------------------------------%
%---------Calculate and Display Results---------%
%-----------------------------------------------%


SL         = zeros(length(y), 1);
SL_overlap = zeros(length(y), 1);

for j = 1:1:length(SL)
    [SL(j) dSL_dt] = Ldep_getSL(T(j), SL_params);
    k = getOverlap(SL(j));
    SL_overlap(j) = k;
end

force1 = SL_overlap(:,1).*(y(:,6).*y(:,8).*(y(:,8)>0) + y(:,7).*y(:,9).*(y(:,9)>0));    % Stu added terms here that zero out tension if the bridges have negative mean distortion

figure(1)
plot(T, force1)

figure(2)
plot(T,SL)











