% this file is published under Creative Commons BY-NC-SA
% 
% Assimilating genome-scale metabolic reconstructions with modelBorgifier
% in preparation
%
% John T. Sauls and Joerg M. Buescher
% BRAIN Aktiengesellschaft
% Microbial Production Technologies Unit
% Quantitative Biology and Sequencing Platform
% Darmstaeter Str. 34-36
% 64673 Zwingenberg, Germany
% www.brain-biotech.de
% jrb@brain-biotech.de
% 
%
function cNum = countC(nowFormula)
%countC finds number of carbons in a formula, if any.
%
% cNum = countC(nowFormula)

% Remove alternate formulas.
nowFormula = regexprep(nowFormula,'/|,*','') ;

% Look for C's
cNum = regexp(nowFormula,'C\d*','match') ;
if isempty(cNum)
    cNum = 0 ;
else
    if length(cNum{1}) == 1 ;
        cNum = 1 ;
    else
        cNum = str2double(cNum{1}(2:end)) ;
    end
end
