function obj = optWeightLin(weight, highscoreVectors, lowscoreVectors)
% this is the objective function for the linear optimization of the weighting of 
% the individual scores. It maximizes the difference between the weighted sum of 
% scores of the correctly assigned reactions and the incorrect matches in the 
% training data set.
%
    obj = 1/(abs(mean(highscoreVectors'*weight) - ...
                 mean(lowscoreVectors'*weight))+1) ;
end
