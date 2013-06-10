function Model = verifyModel(Model)
%verifyModel Ensures that a model is in the correct format to be analyzed
% by any of the scripts in the Tmodel suite.
%   
%INPUTS
% Vmodel    Model from readCbModel any of the readModel functions.
%
%OUTPUTS
% Model     Model with additional fields and correct format.
%
%CALLS
% TmodelFields
% fixNames
% removeDuplicateNames
% makeNamesUnique
% buildRxnEquations
% orderModelFields

%% Declare variables.
nRxns = length(Model.rxns) ;
nMets = length(Model.mets) ;

% Current Field names in the model. 
fieldNames = fieldnames(Model) ;

% Get desired field names. 
fields = TmodelFields ; 
Field.rxn = fields{1} ;
Field.rNum = fields{2} ;
Field.met = fields{3} ;
Field.mNum = fields{4} ; 
Field.all = fields{5} ;
              
%% Pad fields with empty strings or zeros.
% Reaction related cell arrays.
for iField = 1:length(Field.rxn)
   fieldIndex = strcmp(Field.rxn{iField},fieldNames) ;
   fieldIndex = find(fieldIndex) ;
   if isempty(fieldIndex) % Field does not currently exist
       fprintf(['Array .' Field.rxn{iField} ' not in Model. Adding.\n'])
       Model.(Field.rxn{iField}) = cell(nRxns,1) ;
       Model.(Field.rxn{iField})(:) = {''} ;
   else % Field exists, check if the length is correct
       vFieldLength = length(Model.(fieldNames{fieldIndex})) ;
       if vFieldLength ~= nRxns ;
           vFieldLength = vFieldLength + 1 ;
           Model.(fieldNames{fieldIndex})(vFieldLength:nRxns) = {''} ;
       end
   end
end

% Reaction related double arrays.
for iField = 1:length(Field.rNum)
   fieldIndex = strcmp(Field.rNum{iField},fieldNames) ;
   fieldIndex = find(fieldIndex) ;
   if isempty(fieldIndex) % Field does not currently exist
       fprintf(['Array .' Field.rNum{iField} ' not in Model. Adding.\n'])
       Model.(Field.rNum{iField}) = zeros(nRxns,1) ; 
   else % Field exists, check if the length is correct
       vFieldLength = length(Model.(fieldNames{fieldIndex})) ;
       if vFieldLength ~= nRxns ;
           vFieldLength = vFieldLength + 1 ;
           Model.(fieldNames{fieldIndex})(vFieldLength:nRxns) = 0 ;
       end
   end
end

% Metabolite related cell arrays.
for iField = 1:length(Field.met)
   fieldIndex = strcmp(Field.met{iField}, fieldNames) ;
   fieldIndex = find(fieldIndex) ;
   if isempty(fieldIndex) % Field does not currently exist
       fprintf(['Array .' Field.met{iField} ' not in Model. Adding.\n'])
       Model.(Field.met{iField}) = cell(nMets,1) ;
       Model.(Field.met{iField})(:) = {''} ;
   else % Field exists, check if the length is correct
       vFieldLength = length(Model.(fieldNames{fieldIndex})) ;
       if vFieldLength ~= nMets ;
           vFieldLength = vFieldLength + 1 ;
           Model.(fieldNames{fieldIndex})(vFieldLength:nRxns) = {''} ;
       end
   end
end

% Metabolite related double arrays.
for iField = 1:length(Field.mNum)
   fieldIndex = strcmp(Field.mNum{iField}, fieldNames) ;
   fieldIndex = find(fieldIndex) ;
   if isempty(fieldIndex) % Field does not currently exist
       fprintf(['Array .' Field.mNum{iField} ' not in Model. Adding.\n'])
       Model.(Field.mNum{iField}) = zeros(nRxns,1) ; 
   else % Field exists, check if the length is correct
       vFieldLength = length(Model.(fieldNames{fieldIndex})) ;
       if vFieldLength ~= nRxns ;
           vFieldLength = vFieldLength + 1 ;
           Model.(fieldNames{fieldIndex})(vFieldLength:nRxns) = 0 ;
       end
   end
end

% Redeclare field names
fieldNames = fieldnames(Model) ;


%% Ensure reactions are forwards.
fprintf('Making sure reactions are all forwards\n')
% Find reverse reactions based on bounds.
revRxns = find(abs(Model.lb) > Model.ub) ; 

% Do it.
for iRxn = 1:length(revRxns)
    % Reverse bounds.
    [Model.lb(revRxns(iRxn)),Model.ub(revRxns(iRxn))] = ...
        deal(-Model.ub(revRxns(iRxn)), -Model.lb(revRxns(iRxn))) ;
    
    % Change sign of stochiometrix matrix.
    metStoics = find(Model.S(:,revRxns(iRxn))) ;
    for iMet = 1:length(metStoics)
        if Model.S(metStoics(iMet),revRxns(iRxn)) > 0
            Model.S(metStoics(iMet),revRxns(iRxn)) = ...
                0 - Model.S(metStoics(iMet),revRxns(iRxn)) ;
        else
            Model.S(metStoics(iMet),revRxns(iRxn)) = ...
                abs(Model.S(metStoics(iMet),revRxns(iRxn))) ;
        end
    end
end

%% Remove hyphens and obtuse characters from names.
fprintf('Fixing names of metabolites and reaction\n')
Model.mets = fixNames(Model.mets) ;
Model.metNames = fixNames(Model.metNames) ;
Model.metNames = removeDuplicateNames(Model.metNames) ; 
Model.rxns = fixNames(Model.rxns) ;
Model.rxnNames = removeDuplicateNames(Model.rxnNames) ; 

%% Make sure all reaction and metabolite names are unique.
fprintf('Checking if reaction IDs (.rxns) are unique.\n') ; 
if length(Model.rxns) ~= length(unique(Model.rxns))
    fprintf('ERROR: Not all reactions are unique.\n')
    % Launch name correcting script.
    Model.rxns = makeNamesUnique(Model.rxns,Model.rxnNames) ; 
end
fprintf('Checking if metabolite IDs (.mets) are unique.\n') ;
if length(Model.mets) ~= length(unique(Model.mets))
    fprintf('ERROR: Not all metabolites are unique.\n')
    Model.mets = makeNamesUnique(Model.mets,Model.metNames) ; 
end

%% Make sure all metabolites have a compartment.
noComp = find(cellfun(@isempty,regexp(Model.mets,'\[\w\]$'))) ;

if isempty(noComp)
    fprintf('All metabolites have comparment designation.\n')
else
    % try to get compartment information from metNames
    nameComp = find(~cellfun(@isempty,regexp(Model.metNames,'\[\w\]$'))) ;
    for inc = rowvector(nameComp)
        Model.mets{inc} = [Model.mets{inc} ...
            Model.metNames{inc}(strfind(Model.metNames{inc},'['): ...
            strfind(Model.metNames{inc},']')) ] ;
        Model.metNames{inc} = ...
            Model.metNames{inc}(1:strfind(Model.metNames{inc},'[')-1) ;
    end
    
    noComp = find(cellfun(@isempty,regexp(Model.mets,'\[\w\]$'))) ;

    if isempty(noComp)
        fprintf('All metabolites have comparment designation.\n')
    else    
        fprintf(['Some metabolites have no compartment designation. ' ...
            'Assigning as cytosolic\n'])
        for iMet = 1:length(noComp)
            Model.mets{noComp(iMet)} = [Model.mets{noComp(iMet)} '[c]'] ;
        end
    end
end

%% Give .rxnID or .metID names from .rxns and .mets.
Model.rxnID = Model.rxns ;
Model.metID = Model.mets ;

%% Format chemical formulas
Model.metFormulas = fixChemFormulas(Model.metFormulas) ;

%% Rebuild reaction equations to ensure they use fixed names. 
Model = buildRxnEquations(Model) ; 

%% Ensure all vectors are column vectors.
colFields = {'rxn' 'rNum' 'met' 'mNum'} ;
for iName = 1:length(colFields)
    for iField = 1:length(Field.(colFields{iName})) ;
        Model.(Field.(colFields{iName}){iField}) = ...
            Model.(Field.(colFields{iName}){iField})(:) ;
    end
end
             
%% Make sure there is a model name and if not, ask for one.
if ~isfield(Model,'description')
    needModelName = 1 ; 
else
    fprintf(['Current model name is:\n' Model.description '\n'])
    answer = input('Keep name? (y/n): ','s') ;
    if strcmpi(answer,'y') || strcmpi(answer,'yes')
        fprintf('Keeping name\n')
        needModelName = 0 ;
    else
        needModelName = 1 ; 
    end
end
if needModelName
    prompt = 'Input model name: ' ;
    modelName = input(prompt,'s') ;
    Model.description = modelName ;
end

%% Reorder fields and organize model based on most common mets.
Model = orderModelFields(Model) ; 
Model = organizeModelCool(Model) ; 
