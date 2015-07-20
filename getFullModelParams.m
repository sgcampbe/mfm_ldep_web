%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: getFullModelParams                                        %
%   Date Started: 8/25/2008                                             %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function returns all necessary parameters to
%   perform simulations using the full model. It concatenates structs
%   containing iRU and eRU parameters.
%   See Program Glossary for variable definitions.
%
%   4/14/09 - This version has been updated for the 4 state iRU model.
%
%   7/15/09 - Updated to 5 states in iRU model.
%***********************************************************************%

function param_out = getFullModelParams(numRUs, bparams)

iRUparams = getiRUParams(bparams);

eRUparams = geteRUParams(numRUs, bparams);

% Concatenate iRU and eRU parameters into a single struct
names     = fieldnames(iRUparams);

for i = 1:length(names)
    param_out.(names{i}) = iRUparams.(names{i});
end

names     = fieldnames(eRUparams);

for i = 1:length(names)
    param_out.(names{i}) = eRUparams.(names{i});
end

return