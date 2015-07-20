%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: solveSSiRU                                                %
%   Date Started: 8/26/2008                                             %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function solves the iRU model for a given input Ca
%   concentration.
%   Note that AiRU is the full, multi-D version that has not yet been
%   collapsed and finalized.
%
%   Note that this uses the neighbor-independent values for kon and koff,
%   so it truly is the INDIVIDUAL RU model.  In the coupled iRU/eRU model,
%   the result xiRU has to be altered further, which is done in
%   solveSingleFullSS.m.
%***********************************************************************%

function xiRU = solveSSiRU(AiRU, namesp, Ca_i, xMpr, xMpo, x0)


kbdex  = strcmp('kb',namesp);                   % Index to kb matrix within AiRU
hfdex  = strcmp('hf',namesp);                   % Index to hf matrix within AiRU
gxbdex = strcmp('gxb',namesp);                  % Index to gxb matrix within AiRU
Akb    = AiRU(:,:,kbdex);                       % Extract matrices from AiRU
Ahf    = AiRU(:,:,hfdex);                       % 
Agxb   = AiRU(:,:,gxbdex);                      % 
summap = ~(kbdex + gxbdex);                     % Logical indexing structure to sum all other matrices into final AiRU
% summap = ~(kbdex+hfdex+gxbdex);                 % Logical indexing structure to sum all other matrices into final AiRU
AiRU   = sum(AiRU(:,:,summap),3);               % Sum all kinetic matrices save Akb 

AiRU   = updateAiRU(AiRU, ...                   % Update kinetic matrix
                    Akb,  ...
                    Ca_i,  ...
                    Agxb,  ...
                    xMpr,  ...
                    xMpo,  ...
                    x0);


% AiRU   = updateAiRU(AiRU, ...                   % Update kinetic matrix
%                     Akb,  ...
%                     Ca_i, ...
%                     Ahf, ...
%                     Agxb, ...
%                     xPr, ...
%                     xPo, ...
%                     x0);
                
xiRU   = solveEigSS(AiRU);                      % Solve steady-state

return