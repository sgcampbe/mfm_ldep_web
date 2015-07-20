%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: getiRUParams                                              %
%   Date Started: 8/9/2008                                              %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function sets up the individual RU model kinetic
%   coefficient matrix for the simple 4 state iRU model.  The approach of
%   building perfectly generic matrices is kind of overkill, but in order
%   to preserve the same basic structure, this approach has been retained.
%
%   7/14/09 - Major updates to add 3 state XB model with strain-dependent
%   kinetic parameters.
%***********************************************************************%

function param_out = getiRUParams(p)

AiRU  = zeros(4, 4, 6);       % Initialize output

% Set up connectivity graphs or maps linking iRU model states
% The order of setup should correspond with that defined in getParams.m
map.kb    =  [ 1,  2];
         
map.ku    =  [ 2,  1];

map.kon   =  [ 2,  3];
             
map.koff  =  [ 3,  2]; 
            
map.f     =  [ 3,  4]; 
            
map.g     =  [ 4,  3];

map.hf    =  [ 4,  5];

map.hb    =  [ 5,  4];

map.gxb   =  [ 5,  3];

names  = fieldnames(map);    % Retrieve field names for maps
namesp = fieldnames(p);      % Retrieve field names for parameter struct
            
% Check to assure map order corresponds to order defined in p:
trues = sum(strcmp(names, namesp(1:length(names))));
if trues < length(names)
    error('Parameter order in this function does not match that of getParams!');
end

%-----------------------%
% Create AiRU 3D matrix %
%-----------------------%

% Create kinetic matrix using maps and corresponding parameter values
for n = 1:length(names)
    I = map.(names{n})(2);
    J = map.(names{n})(1);
    for i = 1:length(I)
        AiRU(I(i),J(i), n) = p.(names{n});
        AiRU(J(i),J(i), n) = -p.(names{n});
    end
end

% Set up output
param_out.AiRU      = AiRU;
param_out.namesp    = namesp(1:9);         % NOTE: Only names of AiRU 'layers' are retained

return


    

