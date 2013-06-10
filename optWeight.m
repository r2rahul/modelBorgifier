function obj = optWeight(weight, highscoreVectors, lowscoreVectors)
    obj = 1/(abs(mean(highscoreVectors'*weight) - ...
                 mean(lowscoreVectors'*weight))+1) ;
end