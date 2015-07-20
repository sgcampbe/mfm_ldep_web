cls

colors = struct('wt', 'k', ...
                'het', 'r', ...
                'hetcomp', 'r--');


% Preliminary declarations...
fit_tau_flag = false;   % Set to true if any experimental data changes...
ca_data_dir  = fullfile('PKD2_simulations', 'Ca_data');   % Directory containing Ca data (relative path)
numRUs       = 13;
t_end        = 1000;        % [ms]
task_params.L_code   = 1;           % Isometric total muscle length - just a formality since series stiffness is zero
task_params.L_min    = 4;           % Doesn't really matter - just needs to be longer than SL_slack, probably


% Add model directory to the path
cd('..')
path(fullfile('.', 'PKD2_simulations'), path);

%-------------------------------------%
% Enter Average Data from Experiments %
%-------------------------------------%

rmag    = struct('wt',  0.3, ...
                 'het', 0.37);

tau     = struct('wt',  170.2, ...
                 'het', 185.1);

Ca_ttp  = struct('wt',  54.9, ...
                 'het', 53.0);

RT50    = struct('wt',  66.9, ...
                 'het', 66.0);

SL_peak = struct('wt',  3.8, ...
                 'het', 4.0);

SL_ttp  = struct('wt',  103.3, ...
                 'het', 101.8);
             
% Assume absolute Ca values:
Ca_min     = 0.1;   % µm
Ca_max_ref = 1.0;   % µm, associated with WT tyrodes for reference

stim       = 50;  % seconds, stimulus time offset

%-------------------------------------------%
% Make Synthetic Ca transients - WT and Het %
%-------------------------------------------%

if fit_tau_flag
    strains = {'wt', 'het'};
    for i = 1:length(strains)
        strain = strains{i};

        % Estimate absolute Ca max value based on Rmag
        Ca_max = Ca_min + rmag.(strain) * ((Ca_max_ref - Ca_min) / rmag.wt);

        [t_data, Ca_data, tau_1, tau_2] = fitSyntheticCaT(Ca_ttp.(strain), tau.(strain), Ca_min, Ca_max, stim);

        % Save the file containing Ca_data and t_data
        ca_file_str = sprintf('CaT_%s', strain);
        save(fullfile(ca_data_dir, ca_file_str), 't_data', 'Ca_data', 'tau_1', 'tau_2')
    end
end

%-------------------------------------%
% Define Model Parameters for Fitting %
%-------------------------------------%



%----------------------------------------------%
% Run Three Simulations: WT, HET, and HET-pTnI %
%----------------------------------------------%

sims = {'wt', 'het', 'hetcomp'};
ca_name = {'wt', 'het', 'het'};


% Set up plots
axh = multiSubplot(2,4);

% Set up parameter values for the different simulations

basickey  = {};
basicvals = {}; 

% Add other mutual parameter values
argkey     = {'p_f', 'eps_RLCP', 'k_s', 'Temp', 'fapp', 'gxb', 'hf', 'hb', 'xbmodsp', 'kon', 'gammaB', 'kxb', 'ku'};
argvals.wt = {1.632,   2.462,      0,   273+30,  0.10,   1.1,  2.0,   0.4     1.0,     750,    300      40   0.625};

ku_tniphos = 0.775; % [ms] Ca dissociation rate under elevated TnI phosphorylation

argvals.het = argvals.wt;
argvals.hetcomp = argvals.wt;

argvals.hetcomp{end} = ku_tniphos;


% Start Sim Loop
for i = 1:length(sims)
    
    sim = sims{i};
    
    %----------------------------%
    %  Set Ca Transient driver   %
    %----------------------------%

    Ca_params = {0, 'none', 'none', fullfile(ca_data_dir, ['CaT_' ca_name{i}])};
    L_params  = {task_params.L_code task_params.L_min task_params.L_min 0 0};   % Codes for unloaded cell shortening


    %-----------------------------------%
    %   Get Auto-generated model data   %
    %-----------------------------------%

    bparams = getParams(argkey, argvals.(sim));         % Get basic model params
    params  = getFullModelParams(numRUs, bparams);  % Get all auto generated structures


    %***************************%
    %   Run Twitch Simulation   %
    %***************************%

    tic
    rawdata = twitchResponse(params, t_end, Ca_params, L_params);
    toc


    %--------------------------%
    %   Data Post-processing   %
    %--------------------------%
    force = calcForce(rawdata.X, params, rawdata.T, L_params);
    t.(sim)     = rawdata.T;
    SL.(sim)    = rawdata.SL;
    Ca.(sim)    = rawdata.Ca;
    
    % Add plot lines
    
    axes(axh(1,1))
    hold on
    plot(t.(sim), Ca.(sim), colors.(sim))

    axes(axh(2,1))
    hold on
    plot(t.(sim), SL.(sim), colors.(sim))
    
    
    props.(sim) = getShorteningProps(t.(sim), SL.(sim));

end

%------------------%
% Plot the Results %
%------------------%

axes(axh(1,2))
xlabel(axh(1,2), 'WT Peak SL Short.')
plot([1 2], [SL_peak.wt props.wt.SH_peak], 'b^')

axes(axh(1,3))
xlabel(axh(1,3), 'WT SL_ttp')
plot([1 2], [SL_ttp.wt props.wt.t_peak], 'b^')

axes(axh(1,4))
xlabel(axh(1,4), 'WT RT50')
plot([1 2], [RT50.wt props.wt.t_pt50], 'b^')

axes(axh(2,2))
xlabel(axh(2,2), 'HET Peak SL Short.')
plot([1 2 3], [SL_peak.het props.het.SH_peak props.hetcomp.SH_peak], 'r^')

axes(axh(2,3))
xlabel(axh(2,3), 'HET SL_ttp')
plot([1 2 3], [SL_ttp.het props.het.t_peak props.hetcomp.t_peak], 'r^')

axes(axh(2,4))
xlabel(axh(2,4), 'HET RT50')
plot([1 2 3], [RT50.het props.het.t_pt50 props.hetcomp.t_pt50], 'r^')



%--------------%
% Save Results %
%--------------%

