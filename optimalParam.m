function [bestCost bestGamma] = optimalParam(grouplabels,data)
%optimalParam determines the best cost and gamma parameters for svmtrain. 
%
% [bestCost bestGamma] = optimalParam(grouplabels,data)
%
%INPUTS
% grouplabels   Array of hits and misses (ones and zeros), corresponding to
%               the data in data.
% data          Array of data to train against. 
%
%OUTPUTS
% bestCost      Cost parameter which leads to the heigh
% bestGamma
%
%CALLED BY
% optimalScores

stepsize = 10 ;
costrange = -30:stepsize:30 ;
gammarange = -20:stepsize:40 ;
lastqc = 0 ;
bestqc = 0 ;
improve = true ;
backupcnt = 0 ;
while improve
    improve = false ;
    backupcnt = backupcnt + 1 ;
    clear nowqc ; clear qc ;
    ic = 0 ;
    for cost = costrange
        ic = ic + 1 ;
        ig = 0 ;
        for gamma = gammarange
            ig = ig + 1 ;
            nowqc(ic,ig) = svmtrain(grouplabels, data, ...
                ['-v 3 -c ' num2str(2^cost) ' -g ' num2str(2^gamma) ...
                 ' -q -p .1']) ;
            if nowqc(ic,ig) > bestqc
                bestCost = cost ;
                bestGamma = gamma ;
            end
        end
    end
    
    for ic2 = 2:ic-1
        for ig2 = 2:ig-1
            qc(ic2-1,ig2-1) = mean(mean(nowqc(ic2-1:ic2+1,ig2-1:ig2+1)));
        end
    end
    bestqc2 = max(max(qc)) ;
    bestpos = qc == bestqc2 ;
    bestpos = [round(mean(find(sum(bestpos,2)))) ...
               round(mean(find(sum(bestpos,1)))) ] + 1 ;
           
%     disp([bestCost bestGamma])
    
    % Make search grid finer. 
    stepsize = stepsize * 0.5 ;
    costrange = costrange(bestpos(1)-1):stepsize:costrange(bestpos(1)+1) ;
    gammarange = gammarange(bestpos(2)-1):stepsize:gammarange(bestpos(2)+1) ;
    
    if lastqc < bestqc2 - 0.03
        improve = true ;
        lastqc = bestqc ;
    end
    if backupcnt > 5
        improve = false ;
%         figure
%         surf(gammarange, costrange, nowqc)
%         xlabel('Gamma Range') ; ylabel('Cost Range') ; 
%         zlabel('Accuracy') ; 
    end
end