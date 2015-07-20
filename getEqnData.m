function eqnData = getEqnData(numRUs) 

% Cell array of data structure names
varname = {'alphas', ...
           'm', ...
           'n', ...
           'iPre', ...
           'iPost', ...
           'transTypes'};

eqnData.i_noRUmask = {};    % Initialize to empty cell array
       
for i = 1:length(varname)
    filename = fullfile('EqnData', [varname{i}, num2str(numRUs), '.txt']);
    eval(strcat('eqnData.', varname{i}, '=load(filename);'));
end

return
