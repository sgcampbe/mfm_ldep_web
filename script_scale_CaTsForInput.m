%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Script: script_scaleCaTsForInput.m                                  %
%   Date Started: 3/16/2010                                             %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This script loads experimental records of Ca transients,
%   processes them, scales them (according to estimates from the
%   literature regarding calibrated amplitude), decimates them to
%   reasonable size, and packages them so they're ready to be pasted into
%   the function getCa.m for use in simulations.  This looks at a 2Hz beat
%   and a 4 Hz beat.
%***********************************************************************%



cls


path('/Users/stuart/twitch_dev/MATLAB_scripts', path)   % Get access to fancy functions...

% Make Time vector (we know it's 2 kHz sampling for 0.5 sec)
t = 0:0.5:499.5;

prep = '100302_m01_WT';
    
basedir  = '/Users/stuart/data/PapillaryData_GEN/';
basetag  = 'ISO_2HZ_0MS';
tagsuff  = '_ISO';
diast_autof_flag = false;

prepdir  = [basedir prep '/exported_data'];

%----------------------%
% Import 2Hz transient %
%----------------------%

alldata  = importAllPrepData(prepdir, 0.05, 'end');
[Yratio, Y340, Y380] = calcNormRatio(alldata, basetag, tagsuff, diast_autof_flag);

% Apply some heavy smoothing
Yratio_2hz = sgolayfilt(Yratio,2,41);

%----------------------%
% Import 4Hz transient %
%----------------------%

basetag  = 'ISO_4HZ_0MS';
alldata  = importAllPrepData(prepdir, 0.05, 'end');
[Yratio, Y340, Y380] = calcNormRatio(alldata, basetag, tagsuff, diast_autof_flag);

% Apply some heavy smoothing
Yratio_4hz = sgolayfilt(Yratio,2,41);

figure
plot(Yratio_2hz)
hold on
plot(Yratio_4hz,'r')

%-----------------------------------------------------------%
% Scale transients according to Gao et al, 1998 calibration %
%-----------------------------------------------------------%

% Step 1 - determine scaling based on 2 Hz transient and Gao reported
% systolic and diastolic values
Ca_s = 1.75;    % uM
Ca_d = 0.25;    % uM

Y_s  = max(Yratio_2hz);
Y_d  = mean(Yratio_2hz(1:40));

scale = (Ca_s - Ca_d) / (Y_s - Y_d);

% Step 2 - Scale 2Hz transient and then find necessary offset
Ca_2hz = scale * Yratio_2hz;

offset = mean(Ca_2hz(1:40)) - Ca_d;

% Step 3 - Apply offset to 2hz trace
Ca_2hz = Ca_2hz - offset;

% Step 4 - Scale and offset the 4 Hz trace
Ca_4hz = Yratio_4hz * scale - offset;

% Step 5 - Decimate signals
Ca_2hz = decimate(Ca_2hz, 10);
Ca_4hz = decimate(Ca_4hz, 10);
t      = decimate(t, 10);

% Plot final result
figure
plot(Ca_2hz)
hold on
plot(Ca_4hz, 'r')
