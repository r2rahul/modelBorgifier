function [matchScores,matchIndex,varargout] = findMetMatch(cMet,varargin)
%findMetMatch Finds find possible matches for a given metabolite.
%
% [matchScores,matchIndex[,hit]] = findMetMatch(cMet[,tRxn])
%
%
%INPUTS
% cMet      Metabolite number (relative to CMODEL) metabolite to be
%           compared.
%
%OPTIONAL INPUTS
% tRxn      The reaction from TMODEL that this metabolite's reaction is
%           matched.
%
%GLOBAL INPUTS
% CMODEL
% TMODEL
%
%OUTPUTS
% matchScores  Array of sorted, normalized scores for each met in TMODEL.
% matchIndex   Index of match in TMODEL.mets.
%
%OPTIONAL OUTPUTS
% hit          Flag for if a match was found. 
%
%CALLED BY
% metCompare
% metCompareGUI
%
%DEVELOPMENT NOTES/TO DO
% There is no comparison between CheBIID, PubChemID, or InChIString yet.
% 2013.02.21 Added global functionality


%% Declare variables and scoring structure.
global CMODEL TMODEL 

% Declare matched reaction, if any.
if ~isempty(varargin)
    tRxn = cell2mat(varargin(1)) ;
else
    tRxn = 0 ;
end

% Score values.
ScoreVal = struct('ID', [1,-0.1], ...
    'compartment', [0.01, -5], ... % wrong compartment will be down-weighted again further down
    'name', [1,-0.01], ...
    'formula', [1,-0.0], ...
    'charge', [.1,-.1], ...
    'KEGGID', [1,-0.01], ...
    'SEEDID', [1,-.01], ...
    'ChEBIID', [1,-0.01], ...
    'PubChemID', [1,-0.01], ...
    'InChIString', [1,-0.01]) ;

% Score array.
metScores = zeros(1,length(TMODEL.mets));

% Variable indicate if any significant matching information was found, 
% ie name, formula, ID. If not, this flag indicates a new met should be
% declared. 
hit = 0 ;

% Compartments of every met in T, used later for comparison.
tComp = cell(length(TMODEL.mets),1) ;
for iMet = 1:length(TMODEL.mets)
    tComp{iMet} = TMODEL.mets{iMet}(end-2:end) ;
end


%% Compare cMet to TMODEL
% metsID Match. There is always one ID.
% Match name, removing '[.]' and add regular expression info. 
if ~isempty(regexp(CMODEL.mets{cMet},'\[.\]','once')) ;
    name =  strcat('(\||\<)',CMODEL.mets{cMet}(1:end-3),'(\>|\[|_|\|)') ;
else
    name = strcat('(\||\<)',CMODEL.mets{cMet},'(\>|\[|_|\|)') ;
end
% Match against IDs for all models.
match = regexpi(TMODEL.metID,name) ;
match = ~cellfun(@isempty,match) ;
% If at least one match was found, change flag and give scores.
if ~isempty(find(match,1)) ; 
    hit = 1 ;
    metScores(match) = metScores(match) + ScoreVal.ID(1) ;
end

% Match against compartment for all models.
match = strcmp(tComp,CMODEL.mets{cMet}(end-2:end)) ;
% If at least one match was found, change flag and give scores.
% % % if ~isempty(find(match,1)) ; 
    metScores(match) = metScores(match) + ScoreVal.compartment(1) ;
    metScores(~match) = metScores(~match) + ScoreVal.compartment(2) ;
% % % end

% metNames Match. There can be multiple names.
fullname = CMODEL.metNames{cMet} ;
pipePos = [0 strfind(fullname,'|') length(fullname)+1 ] ;
matchSum = zeros(length(TMODEL.mets),1) ;
for i = 1:length(pipePos) - 1
    name = fullname(pipePos(i)+1:pipePos(i+1)-1) ;
    name = strcat('(\||\<)',name,'(\>|\|)') ;
    match = regexpi(TMODEL.metNames,name) ;
    match = ~cellfun('isempty',match) ;
    matchSum = matchSum + match ;
end
match = logical(matchSum) ;
if ~isempty(find(match,1)) ; 
    hit = 1 ;
    metScores(match) = metScores(match) + ScoreVal.name(1) ;
end

% metFormula Match. There will always be one Formula.
name = strcat('(\||\<)',CMODEL.metFormulas{cMet},'(\>|\|)') ;
match = regexpi(TMODEL.metFormulas,name) ;
match = ~cellfun(@isempty,match) ;
if ~isempty(find(match,1)) ; 
    hit = 1 ;
    metScores(match) = metScores(match) + ScoreVal.formula(1) ;
end

% metCharge match. There will always be a charge.
for i = 1:length(TMODEL.mets)
    if CMODEL.metCharge(cMet) == TMODEL.metCharge(i)
        metScores(i) = metScores(i) + ScoreVal.charge(1) ;
    end
end

% metKEGGID Match. There may multiple KEGGIDs or no KEGGID.
if ~isempty(CMODEL.metKEGGID{cMet})
    fullname = CMODEL.metKEGGID{cMet} ;
    pipePos = [ 0 strfind(fullname,'|') length(fullname)+1 ] ;
    matchSum = zeros(length(TMODEL.mets),1) ;
    for i = 1:length(pipePos) - 1
        name = fullname(pipePos(i)+1:pipePos(i+1)-1) ;
        match = regexpi(TMODEL.metKEGGID,name) ;
        match = ~cellfun('isempty',match) ;
        matchSum = matchSum + match ; 
    end
    match = logical(matchSum) ; 
    if ~isempty(find(match,1)) ;
        hit = 1 ;
        metScores(match) = metScores(match) + ScoreVal.KEGGID(1) ;
    end
end

% metSEEDID Match. There may multiple SEEDIDs or no SEEDID.
if ~isempty(CMODEL.metSEEDID{cMet})
    fullname = CMODEL.metSEEDID{cMet} ;
    pipePos = [ 0 strfind(fullname,'|') length(fullname)+1 ] ;
    matchSum = zeros(length(TMODEL.mets),1) ;
    for i = 1:length(pipePos)-1
        name = fullname(pipePos(i)+1:pipePos(i+1)-1) ;
        match = regexpi(TMODEL.metSEEDID,name) ;
        match = ~cellfun('isempty',match) ;
        matchSum = matchSum + match ; 
    end
    match = logical(matchSum) ; 
    if ~isempty(find(match,1)) ;
        hit = 1 ;
        metScores(match) = metScores(match) + ScoreVal.SEEDID(1) ;
    end
end


%% Normalize scores.
% remove negative scores
metScores(metScores < 0) = 0 ;

if max(metScores) > 0 % Avoid NaN.
    metScores = metScores - min(metScores) ;
    metScores = metScores./max(metScores) ;
else
    metScores = 0.01*ones(size(metScores)) ;
end

%% Take some metabolites out of the running.
% If tRxn is provided, then reduce all scores for noninvolved mets by 90%.
if tRxn ~= 0
    metScores(TMODEL.S(:,tRxn) == 0) ...
        = metScores(TMODEL.S(:,tRxn) == 0) * 0.1 ;
end

% Reduce scores of metabolites that don't have the right compartment.
if ~isempty(regexp(CMODEL.mets{cMet},'\[.\]$','once','match')) ;
    nowCom = regexp(CMODEL.mets{cMet},'\[.\]$','once','match') ;
    % Match against .mets, which should contain the comparment designator.
    metScores(~strcmp(nowCom,tComp)) ...
        = metScores(~strcmp(nowCom,tComp)) * 0.1 ;
end


%% Output
[matchScores,matchIndex] = sort(metScores,'descend') ;

varargout = cell(1) ; 
varargout{1} = hit ;

