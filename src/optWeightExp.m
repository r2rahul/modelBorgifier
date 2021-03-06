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
function obj = optWeightExp(weight, highscoreVectors, lowscoreVectors, wnum, hitnum, missnum)
    obj = 1/(abs(mean((highscoreVectors.^repmat(weight(:,2),1,hitnum ))'*weight(:,1)) - ...
                 mean(( lowscoreVectors.^repmat(weight(:,2),1,missnum))'*weight(:,1)))+1) ;
end