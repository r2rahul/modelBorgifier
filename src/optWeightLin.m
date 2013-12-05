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
function obj = optWeightLin(weight, highscoreVectors, lowscoreVectors)
% this is the objective function for the linear optimization of the weighting of 
% the individual scores. It maximizes the difference between the weighted sum of 
% scores of the correctly assigned reactions and the incorrect matches in the 
% training data set.

    obj = 1/(abs(mean(highscoreVectors'*weight) - ...
                 mean(lowscoreVectors'*weight))+1) ;
end
