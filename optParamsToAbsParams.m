function [abs_argvals abs_argkey opt_vals] = optParamsToAbsParams(opt_params, opt_names, param_file_name, varargin)
%************************************************************************%
%   Length-dependent Myofilament Model                                   %
%   optParamsToAbsParams                                                 %
%   Date Started: 12/29/2011                                             %
%   Author: Stuart Campbell                                              %
%   Description: This function turns scaled parameter values such as those
%   used by an optimization routine into absolute parameters based on range
%   definitions from the text file specified by param_file_name.  opt_names
%   is a cell array of strings specifying the model variables contained in
%   opt_params.
%   If opt_params and opt_names do not have the same length or if either is
%   empty, then the parameter values as defined by the param file are used
%   instead.  Any parameter entries in the param file that are not in
%   opt_names are loaded as they appear in the param file, essentially
%   acting as default (non-fitted) values.
%
%
%   12/30/2011 - Adding the ability to optionally save a new param file
%   that updates opt_vals according to the passed in opt_params. This is
%   useful for storing optimization results.  Just supply an additional
%   arguement containing a destination filename for the new file.  It will
%   go into varargin.
%************************************************************************%


%--------------------%
% Read in param file %
%--------------------%

paramstr = file2str(param_file_name);
lines    = split(paramstr, '\n');


%-------------------------------%
% Extract values from each line %
%-------------------------------%

numparams = length(lines);

% Initialize variables
abs_argkey   = cell(1,numparams);
opt_vals     = zeros(1,numparams);  % Scaled parameter values from file
abs_min_vals = zeros(1,numparams);  % Minimum range for absolute values (from file)
abs_max_vals = zeros(1,numparams);  % Maximum range for absolute values (from file)

for i = 1:numparams
    parts = split(lines{i});    % Split by whitespace
    abs_argkey{i} = parts{1};
    abs_min_vals(i) = str2double(parts{2});
    abs_max_vals(i) = str2double(parts{3});
    opt_vals(i)     = str2double(parts{4});
end


%--------------------------------------------------%
% Overwrite opt_vals from any passed-in opt_params %
%--------------------------------------------------%

if (length(opt_params) == length(opt_names)) && ~(isempty(opt_params) || isempty(opt_names))
    numreplacements = length(opt_params);
    for i = 1:numreplacements
        repdex = namesToListIndices(abs_argkey, opt_names{i});
        if isempty(repdex)
            error('opt_name %s does not match any defined parameters in %s', opt_names{i}, param_file_name)
        else
            opt_vals(repdex) = opt_params(i);   % Replace with fitting value
        end
    end
end
        


%--------------------------------%
% Convert opt values to absolute %
%--------------------------------%

% This uses 'triangle' periodicity to force opt values to specify an
% absolute value within the min/max range.  In other words, -0.1 for an opt
% param will be converted to 0.1 and so on.

opt_vals = abs(opt_vals);           % Function is mirrored about zero, so just use abs vals
base     = floor(opt_vals);         % Isolate whole number parts
odds     = boolean(mod(base, 2));   % Identify odd numbers; these need to be 'flipped' to make neg. slope of triangle...
opt_vals = mod(opt_vals, 1);        % reduce to non-whole number parts

opt_vals(odds) = 1 - opt_vals(odds);    % Subtract fractional parts of odd-based opt_vals from one to make triangle

abs_argvals = abs_min_vals + (abs_max_vals - abs_min_vals) .* opt_vals; % Interpolate on range according to opt_vals to find final abs param values

abs_argvals = mat2cell(abs_argvals, 1, ones(1,length(abs_argvals)));    % Convert to cell array of numbers


%----------------------------------------%
% Optionally Save Updated Parameter File %
%----------------------------------------%
if nargin > 3
    template = '%s\t%f\t%f\t%f\n';
    outputfname = varargin{1};
    fid = fopen(outputfname, 'w');
    for i = 1:numparams
        fprintf(fid, template, abs_argkey{i}, abs_min_vals(i), abs_max_vals(i), opt_vals(i));
    end
    fclose(fid);
end


return

