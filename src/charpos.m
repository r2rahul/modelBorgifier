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
function letterpos = charpos(instr)
% finds the letter characters (as opposed to numeric characters) in a string
% and returns a logical array

letterpos = false(size(instr)) ;
for i = 1:length(instr)
    letterpos(i) = isnan(str2double(instr(i))) ;
end
