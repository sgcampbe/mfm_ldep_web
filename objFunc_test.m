function [err comp_vals] = objFunc_F344NIA25C(t, SL, force, task_params)
%-------------------------------------------------------------------------%
%   Objective function for evaluating model fits to F344NIA myocyte data  %
%   Author: Stuart Campbell
%   Date Started: 12/29/2011
%   
%   Description: Objective functions are used in automated fitting
%   processes to compute the error between properties of a twitch 
%   simulation and those of real data.  The simulation data are provided as
%   full timecourses in the vectors t, SL, and force.  task_params is a
%   catch-all struct that contains information such as the file name
%   containing the data to be fitted, etc.
%
%   This function uses calcScaledError.m to compute err, and also packs
%   important info about the comparison into comp_vals.  Lastly, this
%   function plots some output to give the user an idea of how the fit is
%   going.
%-------------------------------------------------------------------------%


%--------------------------------%
% Load Measured Shortening Props %
%--------------------------------%

load(task_params.data_file_name)    % Must load SH_props (struct), stim, and stimtype


%--------------------------------%
% Get Simulated Shortening Props %
%--------------------------------%

[SH_props_sim, ~, markers] = getShorteningProps(t, SL, 'stim', stim, 'stimtype', stimtype, 'sgfilt', 0);


%--------------------------%
% Build Comparison Vectors %
%--------------------------%

propnames = fieldnames(SH_props);   % All properties participating in the fit
numprops  = length(propnames);

guess     = zeros(1,numprops);      % Dimension comparision vectors
actual    = guess;

for i = 1:numprops
    prop = propnames{i};
    guess(i)  = SH_props_sim.(prop);
    actual(i) = SH_props.(prop);
end

%---------------------%
% Calculate Fit Error %
%---------------------%

err = calcScaledError(guess, actual);

err_guess = (actual - guess) ./ actual; % Compute for display purposes.  Scaled, without changes to sign.

%----------------------%
% Pack Comparison Data %
%----------------------%

comp_vals.propnames = propnames;
comp_vals.guess     = guess;
comp_vals.actual    = actual;
comp_vals.t         = t;
comp_vals.SL        = SL;


%---------------%
% Visualize Fit %
%---------------%

f = figure(99);
set(f, 'Position', [396   660   345   438]);

ax(1) = subplot(2,1,1);
cla(ax(1))      % Make sure axes are clear
plot(t, SL)     % Plot sim trace
hold on
markers.draw()  % Add property markers

ax(2) = subplot(2,1,2);
plot([1 numprops], [1 1])           % Unity line
hold on
plot(1:numprops, err_guess, 'bo')   % Scaled errors for each prop
set(ax(2), 'XTick', 1:numprops)
set(ax(2), 'XTickLabel', propnames)
rotateXLabels(ax(2), 45)

return
