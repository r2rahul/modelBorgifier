function [rxnList,metList,Stats] = reactionCompare(CmodelIn,TmodelIn, ...
                                                   scoreIn,varargin) 
%reactionCompare is the function front end for autoMatchReactions and 
% the reactionCompareGUI. Loads CMODEL, TMODEL, Score and ScoreTotal from
% globals if they are not provided. If rxnList is not provided then it is
% generated with autoMatchReactions.
%
%[rxnList metList Stats] = reactionCompare(CmodelIn,TmodelIn,scoreIn, ...
%                                          [rxnList,metList,Stats])
%
%INPUTS
% CmodelIn
% TmodelIn
% scoreIn
%
%OPTIONAL INPUTS
% rxnList    Array pairs reactions in CMODEL with matches from TMODEL or
%            declares them as new.
% metList    Array pairs metabolites in CMODEL with matches from TMODEL,
%            new metabolites are given their new met number in TMODEL.
% Stats      Stats array that contains weighting information from previous
%            scoring work. 
%
%OUTPUTS
% rxnList
% metList
% Stats
%
%CALLS
% optimalScores
% autoMatchReactions
% reactionCompareGUI

%% Declare variables
% Make the inputs variables so they can be easily accessed by downstream
% scripts. None of these variables are edited by any of these scripts. 
global CMODEL TMODEL SCORE
CMODEL = CmodelIn ;
TMODEL = TmodelIn ; 
SCORE = scoreIn ;

% Need scoreTotal now to create rxnList if Stats is not provided.
if nargin <= 5
%     Stats = optimalScores(CMODEL,TMODEL,SCORE) ;
    Stats = optimalScores ;
else
    Stats = varargin{3} ;
end

% Was metList supplied?
if nargin >= 5
    metList = varargin{2} ;
else
    metList = zeros(length(CMODEL.mets),1) ;
end

% How 'bout rxnList?
if nargin > 3
    rxnList = varargin{1} ;
else
    rxnList = ones(length(CMODEL.rxns),1)*-1 ;
    [rxnList,metList] = autoMatchReactions(Stats.scoreTotal, ...
                                           rxnList,metList,0.99,0.1,0.01);
end

%% Manual reaction comparison.

% Create information structure to pass to reactionCompareGUI
InfoBall.rxnList = rxnList ;
InfoBall.metList = metList ; 
InfoBall.CmodelName = CMODEL.description ;
InfoBall.Stats = Stats ;

% Launch GUI.
[rxnList,metList,Stats] = reactionCompareGUI(InfoBall) ; 


