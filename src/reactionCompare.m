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
global CMODEL TMODEL SCORE METNAMEMATCH
CMODEL = CmodelIn ;
TMODEL = TmodelIn ; 
SCORE = scoreIn ;

if isempty(SCORE)
    fprintf('Cannot proceed with reactionCompare becuase SCORE is empty.\n')
    return
end
if isempty(TMODEL)
    fprintf('Cannot proceed with reactionCompare becuase TMODEL is empty.\n')
    return
end
if isempty(CMODEL)
    fprintf('Cannot proceed with reactionCompare becuase CMODEL is empty.\n')
    return
end

% How 'bout rxnList?
if nargin >= 4
    if ~isempty(varargin{1})
        rxnList = varargin{1} ;
    else
        rxnList = ones(length(CMODEL.rxns),1)*-1 ;
    end
else
    rxnList = ones(length(CMODEL.rxns),1)*-1 ;
end

% Was metList supplied?
if nargin >= 5
    if ~isempty(varargin{2})
        metList = varargin{2} ;
    else
        metList = zeros(length(CMODEL.mets),1) ;
    end
else
    metList = zeros(length(CMODEL.mets),1) ;
end

% Need scoreTotal now to create rxnList if Stats is not provided.
%     Stats = optimalScores(CMODEL,TMODEL,SCORE) ;
if nargin >= 6
    if ~isempty(varargin{3})
        Stats = varargin{3} ;
    else
        Stats = optimalScores ;
    end
else
    Stats = optimalScores ;
end

% How about metNameMatch?
if nargin >= 7
    METNAMEMATCH = varargin{4} ;
end
    

%% Manual reaction comparison.

% Create information structure to pass to reactionCompareGUI
InfoBall.rxnList = rxnList ;
InfoBall.metList = metList ; 
InfoBall.CmodelName = CMODEL.description ;
InfoBall.Stats = Stats ;
InfoBall.S = CMODEL.S ;

% Launch GUI.
[rxnList,metList,Stats] = reactionCompareGUI(InfoBall) ; 


