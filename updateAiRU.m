%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: updateAiRU                                                %
%   Date Started: 8/7/2008                                              %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: Updates the kinetic coefficient matrix AiRU according
%   to current values of both Ca_i and the conditional probability psi.
%   This function may be used to update only Akb to current Ca_i by simply
%   passing in zero for Aku and psi (as is done in eigenvalue solves).
%   See Program Glossary for variable definitions.
%
%   4/13/09 - Updated for 4 state iRU model
%
%   1/2/2012 - Updated to implement Rice/Razumova strain dependence.
%
%   1/17/2012 - Added scale factor for distortions
%***********************************************************************%

function AiRU = updateAiRU(AiRU,  Akb, Ca_i, Agxb, Ahf, xMpr, xMpo, varargin)

if nargin > 7           % Must be trying to pass in strain dependent parameters
    bparams = varargin{1}; 
    hfmdc   = bparams.hfmdc;
    sigmap  = bparams.sigmap;
    sigman  = bparams.sigman;
    x0      = bparams.x0;   
else                    % Use some defaults
    hfmdc   = 0;    % These values mean NO strain dependence          
    sigmap  = 0;          
    sigman  = 0;          
    x0      = 0.007;        
end

%---------------------------------%
% Apply Distortion Scaling Factor %
%---------------------------------%

hfmd  = exp(-sign(xMpr)*hfmdc*((xMpr/x0)^2));               % Eqn 22 from Rice.  Slows forward rate when pre-rotation is positive, speeds it up when negative
sigma = (xMpo < x0) * sigman + (xMpo <= x0) * sigmap;       % Update sigma depending on mean distortion (positive or negative)
gxbmd = exp(sigma*((xMpo - x0)^2));                         % Eqn 13 from Razumova paper, makes the rate of cross bridge detachment strain dependent


AiRU = AiRU + Akb * Ca_i + Ahf * hfmd + Agxb * gxbmd;

return