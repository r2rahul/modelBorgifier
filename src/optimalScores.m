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
function Stats = optimalScores(varargin) 
%optimalScore takes the 3D SCORE matrix produced from compareCbModels and
%   weights the scores based on a given training set using either SVM or
%   linear optimization. Colapses SCORE into a cRxnN x tRxnN matrix and
%   also produced the structure Stats, which include the best matches for 
%   each reaction in CMODEL and the corresponding index in TMODEL. rxnList
%   is the training set. If rxnList is not provided, then optimalScores
%   just colapses the 3D SCORE matrix to scoreTotal.
%
% Stats = optimalScore([rxnList, optimizer, optParamFlag])
% Stats = optimalScore([rxnList, 'svm', cost, gamma])
% 
%
%OPTIONAL INPUTS
% rxnList
% optimizer     Either 'svm', 'linear', or 'exp'; the latter two are custom
%               functions contained within the functions optWeightLin and
%               optWeightExp.
% optParamFlag  Flag that initiates paramter optimizer for SVM, calculating
%               cost and gamma.
% cost, gamma   Parameters for svm. 
%
%GLOBAL INPUTS
% CMODEL
% TMODEL
% SCORE
%
%OUTPUTS
% Stats             Structure containing:
%   bestMatch
%   bestMatchIndex
%   weightArray
%   scoreTotal      Readjusted scoreTotal based on training set. 
%
%CALLS
% optWeightLin
% opnWeightExp
% svmtrain
% optimalParam
%
%CALLED BY
% reactionCompare
% autoMatchReactions


%% Declare variables
global SCORE

if nargin > 0
    rxnList = varargin{1} ;
    optimizer = varargin{2} ;
    % Must either delcare to optimize params.
    if nargin == 3
        optParamFlag = 1 ;
    % Or provide cost and gamme, even if optimizer is optFun.
    elseif nargin == 4
        optParamFlag = 0 ;
        cost = varargin{3} ;
        gamma = varargin{4} ;
    end
end

%% Compute SCORE total and stats.
[scoreTotal,Stats] = colapseScore(SCORE) ;

%% Optimize SCORE weights.
if exist('rxnList','var')
    % Determine hits and positions in SCORE.
    fprintf('Constructing training set data.\n')
    hitPos = find(rxnList > 0) ;
    hitPos(:,2) = rxnList(hitPos) ;
    hitVec = zeros(size(SCORE,3),length(hitPos)) ;
    
    % Determine misses. Use the next best SCORE after a hit and the best
    % SCORE from a declared new reaction. 
    missPos = find(rxnList == 0) ;
    missVec = zeros(size(SCORE,3), length(hitPos) + length(missPos)) ;
                 
    % Sorted scoreTotal matrix to find misses. 
    [~, sortI] = sort(scoreTotal,2,'descend') ;
    
    % Find the scores for the hits, and mark next best scores as misses.
    for i = 1:length(hitPos)
        hitVec(:,i) = SCORE(hitPos(i,1),hitPos(i,2),:) ;
        % The miss from top 2 matches. Assumes hit is also in the top 2.
        missSet = setdiff(sortI(hitPos(i,1),1:2),hitPos(i,2)) ;
        for j = 1:1 
            nowLow = SCORE(hitPos(i,1),missSet(j),:) ;
            nowLow = reshape(nowLow,size(nowLow,3),size(nowLow,1)) ;
            missVec(:,i*j) = nowLow ;
        end
    end
    
    % Now add to that the reactions declared as new from the training set.
    for i = 1:length(missPos)
        missSet = sortI(missPos(i),1:1) ;
        for j = 1:1
            nowLow = SCORE(missPos(i),missSet(j),:) ;
            nowLow = reshape(nowLow,size(nowLow,3),size(nowLow,1)) ;
            missVec(:,length(hitPos)+i*j) = nowLow ;
        end
    end
    
    % Choose weighting function.
    fprintf('Optimizing.\n')
    if strcmp(optimizer,'svm')
        % Find optimal best cost and gamma parameters. 
        if optParamFlag
            [cost, gamma] = optimalParam([ones(1,size(hitVec,2)) ...
                                            zeros(1,size(missVec,2))]', ...
                                            [hitVec missVec]') ;      
        end
        disp([cost, gamma])
        
        % Now run svm with determined parameters. 
        fprintf('SVMing.\n')
        svmmodel = svmtrain([ones(1,size(hitVec,2)) ...
                             zeros(1,size(missVec,2))]', ...
                            [hitVec missVec]', ...
                            ['-c ' num2str(2^cost) ...
                             ' -g ' num2str(2^gamma) ...
                             ' -q -p .1']) ;
        weights = svmmodel.SVs' * svmmodel.sv_coef ;
        
    elseif strcmp(optimizer,'RF')
        traindata = [hitVec missVec]' ;
        trainlabel = [true(1,size(hitVec,2)) false(1,size(missVec,2))]' ;
        trunk = zeros( size(traindata,1), size(traindata,2), size(traindata,2)) ;
        % construct square training data
        for it = 1:size(trunk,1)
            trunk(it,:,:) = traindata(it,:)' * traindata(it,:) ;
        end
        % calculate discriminating quality of products of single score
        % dimensions
        qual = zeros(size(traindata,2)) ;
        qualsign = zeros(size(traindata,2)) ;
        for it = 1:size(qual,1) ;
            for it2 = it:size(qual,1) ;
                [~,qual(it,it2)] = ttest2(trunk(trainlabel,it,it2), trunk(~trainlabel,it,it2)) ;
                qualsign(it,it2) = (mean(trunk(trainlabel,it,it2)) < mean(trunk(~trainlabel,it,it2)) ) ;
            end
        end
        % discard useless values
        qual(isnan(qual)) = 1 ;
        qual(qual > 0.1) = 1 ;
        % convert p-value to weight
        weight2 = -log10(qual) ;
        weight2(isinf(weight2)) = 0 ;
        weight2(qualsign == 1) = -weight2(qualsign == 1) ;
        
    elseif strcmp(optimizer,'linear')
        fprintf('Using optWeightLin function.\n')
        weights = fminunc(@(weight)optWeightLin(weight,hitVec,missVec), ...
                  ones(size(hitVec,1),1)) ;
              
    elseif strcmp(optimizer,'exp')
        fprintf('Using optWeightExp function.\n')
        [wnum, hitnum] = size(hitVec) ;
        missnum = size(missVec,2) ;
        signs = sign(mean([hitVec missVec],2)) ;
        weights = fmincon(@(weight)optWeightExp(weight,hitVec,missVec, ...
                                                wnum,hitnum,missnum), ...
                  ones(size(hitVec,1),2), ...
                  [], [], [], [],...
                  [zeros(wnum,1) repmat(0.3,wnum,1)], ...
                  [repmat(1000,wnum,1) repmat(3,wnum,1)] ) ;
        weights(:,1) = weights(:,1) .* signs ;
    elseif strcmp(optimizer,'none')
        weights = ones(size(hitVec,1),1) ;
    end
    
    % Adjust SCORE by the new weights and recompute normalization.
    scoreWeighted = zeros(size(scoreTotal)) ;
    if strcmp(optimizer,'exp')
        for i = 1:size(weights,1)
            scoreWeighted = scoreWeighted + ...
                (double(abs(SCORE(:,:,i))).^weights(i,2)) * weights(i,1) ;
        end
        [scoreTotal,Stats] = colapseScore(scoreWeighted) ; 
        Stats.weights = weights ;
        
    elseif strcmp(optimizer,'RF')
        for it = 1:size(weight2,1)
            for it2 = 1:size(weight2,2)
                scoreWeighted = scoreWeighted + ...
                    double(SCORE(:,:,it)) .* double(SCORE(:,:,it2)) .* weight2(it,it2) ;
            end
        end
        [scoreTotal,Stats] = colapseScore(scoreWeighted) ; 
        Stats.weight2 = weight2 ;
    else
        % Ensure all weights are positive.
        weights(weights < 0) = 0 ;
        
        for i = 1:length(weights)
            scoreWeighted = scoreWeighted + ...
                double(SCORE(:,:,i)) * weights(i) ;
        end
        [scoreTotal,Stats] = colapseScore(scoreWeighted) ;          
        Stats.weights = weights ;
    end
    
    if nargin > 2 && strcmp(optimizer,'svm')
        Stats.cost = cost ;
        Stats.gamma = gamma ; 
    end
end

% Put scoreTotal in Stats for output.
Stats.scoreTotal = scoreTotal ;


%% Subfunctions
% Colapse matrix and normalize scores.
function [scoreTotal,Stats] = colapseScore(scoreWeighted)
global TMODEL CMODEL

% Compute total scores from subscores. If colapseScore is fed an already
% summed SCORE matrix as an argument, this function is impotent.
scoreTotal = sum(scoreWeighted,3) ;
clear scoreWeighted 

% Find the lowest SCORE that isn't a biomass reaction and set the biomass
% reaction scores to that. This is caught by reactions with many mets.
% Also adjust any reaction with more than 10 reactants or products.
scoreTotalHigh = scoreTotal ; 
manyMetC = find(CMODEL.metNums(:,3) > 10) ; 
manyMetC = [manyMetC; find(CMODEL.metNums(:,5) > 10)] ;
scoreTotalHigh(manyMetC,:) = 1000 ;
clear manyMetC
manyMetT = find(TMODEL.metNums(:,3) > 10) ; 
manyMetT = [manyMetT; find(TMODEL.metNums(:,5) > 10)] ;
scoreTotalHigh(:,manyMetT) = 1000 ;
clear manyMetT

% Lowest SCORE which isn't from a biomass reaction.
minScore = min(min(scoreTotalHigh)) ;
% Set found reactions to that.
scoreTotal1k = find(scoreTotalHigh == 1000) ;
clear scoreTotalHigh
scoreTotal(scoreTotal1k) = minScore ; 
clear scoreTotal1k

% Normalize scores.
scoreTotal = scoreTotal + abs(min(min(scoreTotal))) ;
scoreTotal = scoreTotal./max(max(scoreTotal)) ;

% Find best matches.
[bestMatch, bestMatchIndex] = max(scoreTotal,[],2) ;
Stats.bestMatch = bestMatch ;
Stats.bestMatchIndex = bestMatchIndex ;
