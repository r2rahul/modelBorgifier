function obj = optWeightExp(weight, highscoreVectors, lowscoreVectors, wnum, hitnum, missnum)
    obj = 1/(abs(mean((highscoreVectors.^repmat(weight(:,2),1,hitnum ))'*weight(:,1)) - ...
                 mean(( lowscoreVectors.^repmat(weight(:,2),1,missnum))'*weight(:,1)))+1) ;
end