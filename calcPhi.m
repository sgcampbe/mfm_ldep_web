%***********************************************************************%
%   Markov model of thin filament activation                            %
%   Function: calcPhi                                                   %
%   Date Started: 8/7/2008                                              %
%   Author: Stuart Campbell                                             %
%                                                                       %
%   Description: This function calculates the conditional probability
%   P{s2>0|s3=0}, or phi.  Can accept column vectors of inputs and in
%   those cases returns column vectors of phi.
%   See Program Glossary for variable definitions.
%
%   3/17/09 - Altered dynamic calculation of phi to correct for
%   permanently-activated RUs, as is the case when mdex > 1.  Similar to
%   the change required in calcPsi.
%
%   4/16/09 - Altered it again to simply use values of B0 and B1 from xiRU
%   - these have been rigorously reformulated so that they can be used
%   directly in calculating phi and are correct even when NEM-S1 is
%   present.
%***********************************************************************%

function phi = calcPhi(xiRU)

phi = xiRU(2) ./ sum(xiRU(1:2));

return
