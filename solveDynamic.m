%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: solveDynamic                                              %
%   Date Started: 8/7/2008                                              %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: Updates the kinetic coefficient matrix AiRU according
%   to current values of both Ca_i and the conditional probability psi.
%   See Program Glossary for variable definitions.
%
%   4/14/09 - Updated to be specifically for 4 state iRU model and for an
%   eRU model that no longer represents XB binding explicitly.
%
%   NOTE: To simulate a constant SL case, the series stiffness constant in
%   bparams (k_s) must be set to the string 'isometric', which signals that
%   the SL to be used for the simulation will be stored in bparams.SL_slack
%   - that's a little random, but is makes things simpler.
%***********************************************************************%

function data = solveDynamic(x_initial, ... % Vector of initial conditions
                             params,    ... % Struct containing equation data
                             t_start,   ... % [ms] - Time of simulation start
                             t_end,     ... % [ms] - duration of simulation
                             Ca_params, ... % Parameters specifying driving Ca
                             L_params,  ... % Parameters specifying dynamic SL
                             bparams)       % Basic params

% Set integrator options
max_step = 1;   % [ms] - Maximum time step allowed by integrator

%-------------------------------------%
% Extract quantities from params struct
names = fieldnames(params);
for i = 1:length(names)
    eval(strcat(names{i},'=','params.',names{i},';'));
end
% Now all param names are defined
%-------------------------------------%

[xiRU xeRU]     = splitX(x_initial);                   % Separate out xeRU and xiRU

kbdex    = strcmp('kb',namesp);                 % Index to kb matrix within AiRU
Akb      = AiRU(:,:,kbdex);                     % Extract Akb from AiRU
%summap   = ~(kbdex)                            % Logical indexing structure to sum all other matrices into final CAiRU
%AiRU     = sum(AiRU(:,:,summap),3)             % Sum all kinetic matrices save Akb into final AiRU

gxbdex   = strcmp('gxb',namesp);                % Did the same as above but for gxb
Agxb     = AiRU(:,:,gxbdex);

hfdex    = strcmp('hf',namesp);                 % Did the same as above but for hf
Ahf      = AiRU(:,:,hfdex);

summap   = ~(kbdex + gxbdex + hfdex);
AiRU     = sum(AiRU(:,:,summap),3);

SL_slack = bparams.SL_slack;
x0       = bparams.x0;

x_initial = [x_initial; SL_slack; 0; x0];       % Adding intial conditions for SL, xMpr, and xMpo - SL entry is just a dummy (unless bparams.k_s = 'isometric'); has to be determined in next step and updated

if ~ischar(bparams.k_s)                         % As long as this isn't a constant SL simulation...
    SL_0     = calcInitialSL(x_initial, bparams, Ldep_getLtot(0, L_params));
    x_initial(end-2) = SL_0;                    % Update value
end

% Determine Ca tranisent timecourse for simulation
t_Ca    = t_start:((t_end - t_start) / 200):t_end;  % Aim for about 200 points so that the linInterp lookups in getDeriv don't get too costly...
Ca_i    = getCa(t_Ca, Ca_params);

t_span  = t_start:t_end;                        % Vector containing desired timepoints of solution
options = odeset('MaxStep',max_step,'Stats','on');  % Set integrator options
[T,X]   = ode23t(@(t,x) getDeriv(t,         ... % Current time
                                 x,         ... % Vector of state variables
                                 AiRU,      ... % Kinetic matrix for iRU model - UNFINISHED
                                 Akb,       ... % k_b component of AiRU
                                 Agxb,      ... % gxb component of AiRU
                          	     Ahf,       ... % hf component of AiRU
                          	     Akon,      ... % k_on component of AeRU
                         	     Akoff,     ... % k_off component of AeRU
                                 alphas,    ... % Fraction of s3 = 1 RU's in each eRU state
                                 Ca_i,      ... % Vector of Ca transient values
                                 t_Ca,      ... % Vector of time points corresponding to Ca_i vector values
                                 L_params,  ... % Length control parameters to specify dynamic changes in length, if any
                                 bparams),  ... % Basic params
                 t_span, x_initial, options);
             
% X = expandiRU(X, alphas);   % Replace dummy values of [B1] and [C] in xiRU with newly solved ones

[xiRU xeRU] = splitX(X);                    % Split x into iRU and eRU parts    

xMpr = xeRU(:,end-1);                         % Trim out distortion vars
xMpo = xeRU(:,end);


data.SL = xeRU(:,end-2);
data.T  = T;
data.Ca = getCa(T, Ca_params);
data.X  = X;


return