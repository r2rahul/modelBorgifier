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
function Model = organizeModelCool(Model)
%organizeModelCool reorders a COBRA model such that the most promiscuous
% mets are listed first, with sister mets grouped together. Reactions are
% then organized so that the S matrix appoximates a rank order matrix, with
% exchange reactions at the end. 
% 
% Model = organizeModelCool(Model)
%
%OUTPUTS
% Model     Model in COBRA format. Also accept a Tmodel.
%
%INPUTS
% Model     Organized model.
% 
%CALLS
% TmodelFields
%
%CALLED BY
% verifyModel


%% Declare variables for met.
nMets = length(Model.mets) ;
metIndex = zeros(nMets,1) ; 
metPos = 1 ; 

% Unique metabolites sans compartment, used to find sisters.
metsNoComp = cell(nMets,1) ;
for iMet = 1:nMets
    metsNoComp{iMet} = Model.mets{iMet}(1:end-3) ;
end

%% Order metabolites
% Sort by most promiscuous metabolites, grouping sisters.
rxnCount = full(sum((Model.S ~= 0),2)) ;

while metPos <= nMets
    % Most promiscuous met.
    [~,maxIndex] = max(rxnCount) ;

    % If the is a tie, choose the first alphabetically.
    maxIndex = maxIndex(1) ;

    % Find sister metabolites.
    nowMetNameNoComp = Model.mets{maxIndex}(1:end-3) ;
    noCompIndex = find(strcmp(nowMetNameNoComp,metsNoComp)) ;

    % Add mets to future index.
    for iMet = 1:length(noCompIndex)
        metIndex(metPos) = noCompIndex(iMet) ; 
        metPos = metPos + 1 ;
        % Take those mets out of the running. 
        rxnCount(noCompIndex(iMet)) = -1 ; 
    end
end

%% Order reactions. 
nRxns = length(Model.rxns) ;
rxnIndex = zeros(nRxns,1) ; 
% Index and number of exchange reactions.
exIndexes = false(nRxns,1) ;
exIndexes(sum(logical(abs(full(Model.S))),1) == 1) = true ;
nExRxns = sum(exIndexes) ;
exRxnPos = nRxns ; 
rxnPos = nRxns - nExRxns ;
% S matix used to find involved reactions.
SearchS = Model.S ;

% Starting from last metabolite in new order, find involved reactions.
for iMet = nMets:-1:1
    involvedRxns = find(SearchS(metIndex(iMet),:)) ;
    for iRxn = 1:length(involvedRxns)
        % It it is an exchange reaction, put it at the end.
        if exIndexes(involvedRxns(iRxn))
            rxnIndex(exRxnPos) = involvedRxns(iRxn) ;
            exRxnPos = exRxnPos - 1 ;
        else % Otherwise put it here. 
            rxnIndex(rxnPos) = involvedRxns(iRxn) ;
            rxnPos = rxnPos - 1 ;
        end
    end
    % Remove reaction from running.
    SearchS(:,involvedRxns) = 0 ;
end
% put remaining reactions on top
rxnIndex(rxnIndex == 0) = setdiff(1:nRxns,rxnIndex) ;

%% Reorder everything.
% Grab fields names.
fields = TmodelFields ; 

% Reaction related cell field names.
rxnFields = fields{1} ; 
% Reaction related double Fields (which are kept model specific).
rNumFields = fields{2} ;
% Metabolite related cell field names. 
metFields = fields{3} ; 
% Metabolite related double field names. 
mNumFields = fields{4} ; 
% guess correct handling of remaining fields
otherfields = setdiff(fieldnames(Model),[rxnFields; rNumFields; metFields; mNumFields ; {'S' ; 'rxnGeneMat' ; 'description'; 'genes'}]) ;
for io = 1:length(otherfields)
    if all(size(Model.(otherfields{io})) == size(Model.rxns))
        if isnumeric(Model.(otherfields{io}))
            rNumFields{end+1} = otherfields{io} ;
        else
            rxnFields{end+1} = otherfields{io} ;
        end
    elseif all(size(Model.(otherfields{io})) == size(Model.mets))
        if isnumeric(Model.(otherfields{io}))
            mNumFields{end+1} = otherfields{io} ;
        else
            metFields{end+1} = otherfields{io} ;
        end
    else
        warning(['removing handle ' otherfields{io} ])
        Model = rmfield(Model, otherfields{io}) ;
    end
end

% Reorder reaction related lists.
for iField = 1:length(rxnFields)
    if isfield(Model,rxnFields{iField})
        Model.(rxnFields{iField}) = Model.(rxnFields{iField})(rxnIndex) ;
    end
end

% Reorder metabolite related lists.
for iField = 1:length(metFields)
    if isfield(Model,metFields{iField})
        Model.(metFields{iField}) = Model.(metFields{iField})(metIndex) ;
    end
end
for iField = 1:length(mNumFields)
    if isfield(Model,mNumFields{iField})
        Model.(mNumFields{iField}) = ...
            Model.(mNumFields{iField})(metIndex) ;
    end
end

% If this is a Tmodel.
if isfield(Model,'Models')
    % Model names.
    moNames = fieldnames(Model.Models) ;

    % Reorder identity arrays.
    for iMo = 1:length(moNames)
        Model.Models.(moNames{iMo}).rxns = ... 
            Model.Models.(moNames{iMo}).rxns(rxnIndex) ;
        Model.Models.(moNames{iMo}).mets = ...
            Model.Models.(moNames{iMo}).mets(metIndex) ;
    end
    
    for iField = 1:length(rNumFields)
        for iMo = 1:length(moNames)
            Model.(rNumFields{iField}).(moNames{iMo}) = ...
                Model.(rNumFields{iField}).(moNames{iMo})(rxnIndex) ;
        end
    end
else
    for iField = 1:length(rNumFields)
        Model.(rNumFields{iField}) = ...
            Model.(rNumFields{iField})(rxnIndex) ;
    end    
end

% Reorder S matrix.
Model.S = Model.S(:,rxnIndex) ;
Model.S = Model.S(metIndex,:) ;

% Reorder RxnGeneMat
if isfield(Model, 'rxnGeneMat') ;
    Model.rxnGeneMat = Model.rxnGeneMat(rxnIndex,:) ;
end

