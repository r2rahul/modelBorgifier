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
function [TmodelC, Cmodel] = addToTmodel(Cmodel,Tmodel,rxnList,metList, varargin)
%addToTModel Adds reactions and info from a comparison model into the 
% master template database, Tmodel.
%
%Tmodel = addToTModel(Cmodel,Tmodel,rxnList,metList)
%
%INPUTS
% Cmodel        Reaction model to draw new reactions from. Must be in the
%               extended COBRA format.
% Tmodel        Template database in which to add information
%               and new reactions.
% rxnList       Array the length of Cmodel.rxns with corresponding matching
%               reaction numbers from Tmodel. If a rxn is unique in Cmodel,
%               it has 0. rxnList is the output of reactionCompare GUI.
% metList       Array the length of Cmodel.mets with corresponding matching
%               met indexes from Tmodel. New mets should already have their
%               correct index in Tmodel. 
% 'NoClean'     Don't call cleanTmodel at the end of this function
%
%OUTPUTS
% Tmodel        The addended master reaction database.
%
%CALLS
% cleanTmodel
%

%% Declare variables.    
global SCORE METNAMEMATCH
% Number of reactions and metabolites being added. 
cRxns = length(rxnList) ;
cMets = length(metList) ; 

% Cmodel name.
cName = Cmodel.description ;

% Current model names
TmodelNames = fieldnames(Tmodel.Models) ;

% Check if model name already exists in Tmodel.
modelCheck = sum(strcmp(cName,TmodelNames)) ;
while modelCheck > 1
    fprintf(['Cmodel name, ' cName ', already exists in Tmodel.\n'])
    prompt = 'Input new model name for Cmodel: ' ;
    cName = input(prompt,'s') ;
    modelCheck = sum(strcmp(cName,TmodelNames)) ;
end

% Convenient versions for adding information.
cNameNew = strcat(cName,':') ;
cNameAdd = strcat('|',cName,':') ;

% Fields which are kept model specific.
sepFieldNames = ({'rev' 'lb' 'ub' 'c'}) ;

% make sure Tmodel has required fields
if ~isfield(Tmodel,'b')
    Tmodel.b = zeros(size(Tmodel.mets)) ;
end

% Reaction related field names, excluding .rxns, .rxnID and .rxnEquation. 
% anything after getMoDesig does not get a model designation when added.

Tfields = fieldnames(Tmodel) ;
rxnFields = false(size(Tfields)) ;
metFields = false(size(Tfields)) ;
for itf = 1:length(Tfields)
    if size(Tmodel.(Tfields{itf}),1) == size(Tmodel.rxns,1) % if size is that of rxns
        if ~ismember(Tfields{itf}, [sepFieldNames, 'rxns', 'rxnID', 'rxnEquation', 'metNums', 'rxnMetNames' ]) % exclude some fields
            rxnFields(itf) = true ;
        end
    elseif size(Tmodel.(Tfields{itf}),1) == size(Tmodel.mets,1) % if size is that of mets
        if ~ismember(Tfields{itf}, [sepFieldNames, 'mets', 'metID', 'S']) % exclude some fields
            metFields(itf) = true ;
        end
    end 
end
rxnFieldNames = Tfields(rxnFields) ;
metFieldNames = Tfields(metFields) ;

for ir = 1:length(rxnFieldNames)
    if ~isfield(Cmodel,rxnFieldNames{ir})
        if isnumeric(Tmodel.(rxnFieldNames{ir}))
            Cmodel.(rxnFieldNames{ir}) = zeros(size(Cmodel.rxns)) ;
        else
            Cmodel.(rxnFieldNames{ir}) = cell(size(Cmodel.rxns)) ;
        end
    end
end
for ir = 1:length(metFieldNames)
    if ~isfield(Cmodel,metFieldNames{ir})
        if isnumeric(Tmodel.(metFieldNames{ir}))
            Cmodel.(metFieldNames{ir}) = zeros(size(Cmodel.mets)) ;
        else
            Cmodel.(metFieldNames{ir}) = cell(size(Cmodel.mets)) ;
        end
    end
end

getMoDesig = 4 ; 

% Add indentity and number arrays for model being added.
Tmodel.Models.(cName).rxns = false(length(Tmodel.rxns),1) ;
Tmodel.Models.(cName).mets = false(length(Tmodel.mets),1) ;
for iField = 1:length(sepFieldNames)
    Tmodel.(sepFieldNames{iField}).(cName) = 0 ;
end


%% Add genes.
if isfield(Cmodel,'genes') 
    Tmodel.Models.(cName).genes = Cmodel.genes ;
end


%% Deal with new reactions (0 in rxnList).
% Number of reactions currently in Tmodel and in the Tmodel to be.
tRxns = length(Tmodel.rxns) ;
tRxnsNow = tRxns ; 

% Change 0's (new reactions) in rxnList to their new index in Tmodel.
for cRxn = 1:length(rxnList)
   if rxnList(cRxn) == 0
      tRxnsNow = tRxnsNow + 1 ;
      rxnList(cRxn) = tRxnsNow ;
   end
end

if tRxnsNow > tRxns
    % Increase size to Tmodel to accomidate new reactions.
    Tmodel.rxns(end+1:tRxnsNow) = {''} ;
    Tmodel.rxnID(end+1:tRxnsNow) = {''} ;
    Tmodel.rxnEquations(end+1:tRxnsNow) = {''} ;
    for iField = 1:length(rxnFieldNames)
        if iscell(Tmodel.(rxnFieldNames{iField}))
            Tmodel.(rxnFieldNames{iField})(end+1:tRxnsNow) = {''} ;
        end
    end

    % Increase length of identity arrays.
    % Current models.
    for iMo = 1:length(TmodelNames)
        Tmodel.Models.(TmodelNames{iMo}).rxns(end+1:tRxnsNow) = false ; 
    end
    % New model.
    Tmodel.Models.(cName).rxns(end+1:tRxnsNow,1) = true ;

    % Increase size of number array for new 
    for iField = 1:length(sepFieldNames)
        Tmodel.(sepFieldNames{iField}).(cName)(end+1:tRxnsNow,1) = 0 ;
        for iMo = 1:length(TmodelNames)
            Tmodel.(sepFieldNames{iField}). ...
                (TmodelNames{iMo})(end+1:tRxnsNow,1) = 0 ;
        end
    end
    
    % Increase size of comments array.
    Tmodel.rxnNotes(end+1:tRxnsNow,1) = {''} ;

    % Increase size of S matrix
    Tmodel.S(:,end+1:tRxnsNow) = 0 ; 
end

%% Deal with new metabolites (already have index for Tmodel).
% Number of current metabolites in Tmodel.
tMets = length(Tmodel.mets) ; 

% if for some funky reason there are still undecleared metabolites in
% metList, auto-declare them as new. (This really should not happen)
disp([strjoin(Cmodel.mets(metList == 0), ', ') ' was not declared and are treated as new.'])
metList(metList == 0) = max(metList) + (1:sum(metList == 0)) ;

% ensure that there are no gaps in the numbering of new mets
numNewMets = sum(metList > tMets) ;
metList(metList > tMets) = tMets + (1:numNewMets) ;

% Number of new metabolites.
tMetsNow = max(metList); 

if tMetsNow > tMets
    % Increase size to Tmodel to accomidate new metabolites.
    Tmodel.mets(end+1:max(metList)) = {''} ;
    Tmodel.metID(end+1:max(metList)) = {''} ;
    Tmodel.b(end+1:max(metList)) = 0 ;
    for iField = 1:length(metFieldNames)
        if isnumeric(Tmodel.(metFieldNames{iField}))
            Tmodel.(metFieldNames{iField})(end+1:max(metList)) = 0 ;
        else
            Tmodel.(metFieldNames{iField})(end+1:max(metList)) = {''} ;
        end
    end

    % Increase size of identity array are by defnition not in old models. 
    % Current models.
    for iMo = 1:length(TmodelNames)
        Tmodel.Models.(TmodelNames{iMo}).mets(end+1:tMetsNow) = false ; 
    end
    Tmodel.Models.(cName).mets(end+1:max(metList),1) = true ;
    
    % Increase size of S matrix.
    Tmodel.S(end+1:max(metList),:) = 0 ;
end

% Look to see if new metabolites have names that aleady exist in
% Tmodel.mets. Prompt the user to change the name, and also change the name
% for all sister metabolites (the ones that are in other compartments). 
% Change the IDs too.
newMets = find(metList > tMets) ;

% mets from C without comp.
nCmets = length(Cmodel.mets) ; 
metsNoComp = cell(nCmets,1) ;
for iMet = 1:nCmets
    metsNoComp{iMet} = Cmodel.mets{iMet}(1:end-3) ;
end

for iMet = 1:length(newMets)
    % Is there a duplicate?
    if ~isempty(find(strcmp(Cmodel.mets{newMets(iMet)},Tmodel.mets),1))
        % Flag to keep asking for a new name as long as it is not unique
        notUnique = 1 ;
        % Ask for user input
        fprintf(['Metabolite ID %s (%s) already exists in Tmodel, ' ...
                 'with name\n %s.\n'], ...
                Cmodel.mets{newMets(iMet)}, ...
                Cmodel.metNames{newMets(iMet)}, ...
                Tmodel.metNames{find(strcmp(Cmodel.mets{newMets(iMet)}, ...
                                            Tmodel.mets),1)})
        while notUnique
            prompt = 'Rename (w/o compartment) as: ' ;
            newName = input(prompt,'s') ;
            % Check to see if name is already taken.
            if isempty(find(strncmp(newName,Tmodel.mets,length(newName)),1))
                notUnique = 0 ;
            else
                fprintf('Name already taken./n') 
            end
        end
        % Change the name of mets.
        % Find sister mets.
        sisMetPos = find(strcmp(Cmodel.mets{newMets(iMet)}(1:end-3), ...
                                metsNoComp)) ;
        % Change the name of the metabolites
        for iSis = 1:length(sisMetPos) ;
            compartment = Cmodel.mets{sisMetPos(iSis)}(end-2:end) ; 
            Cmodel.mets{sisMetPos(iSis)} = [newName compartment] ; 
        end
    end
end

%% Add reactions.
% find protons and water
ignoreMets = find(strncmpi(Tmodel.mets,'h[',2) + ...
    strncmpi(Tmodel.mets,'h+[',3) + ...
    strncmpi(Tmodel.mets,'proton[',7) +...
    strncmpi(Tmodel.mets,'h2o[',4) + ...
    strncmpi(Tmodel.mets,'water[',6) ) ;

for cRxn = 1:cRxns 
   tRxn = rxnList(cRxn) ;   
   
   % Update identity array
   Tmodel.Models.(cName).rxns(tRxn) = true ;
   % For new reactions, add their name from Cmodel and adjust S matrix. 
   if isempty(Tmodel.rxns{tRxn})
       Tmodel.rxns{tRxn} = Cmodel.rxns{cRxn} ;
       % Add reaction ID
       Tmodel.rxnID{tRxn} = strcat(Tmodel.rxnID{tRxn}, cNameNew, ...
                                   Cmodel.rxnID{cRxn}) ;
       
       % Adjust S matrix
       metPos = find(Cmodel.S(:,cRxn)) ;
       metStoics = Cmodel.S(metPos,cRxn) ;
       for iMet = 1:length(metPos)
           Tmodel.S(metList(metPos(iMet)),tRxn) = metStoics(iMet) ;
       end
   else
       Tmodel.rxnID{tRxn} = strcat(Tmodel.rxnID{tRxn},cNameAdd, ...
           Cmodel.rxnID{cRxn}) ;
   end
   
   % Fill in information arrays for all reactions. This loop ensures that
   % only new information is being added. 
   for iField = 1:length(rxnFieldNames)
     % Current information from C and Tmodel. Makes code cleaner.
     Cinfo = Cmodel.(rxnFieldNames{iField}){cRxn} ;
     Tinfo = Tmodel.(rxnFieldNames{iField}){tRxn} ; 
     % If there is information in the Cmodel field, continue to fill.
     if ~isempty(Cinfo)
       % If there is information already in Tmodel, add to it.
       if ~isempty(Tinfo)
         % Some get model designations.
         if ismember(rxnFieldNames{iField}, {'subSystems', 'rxnECNumbers', 'rxnReferences', 'rxnNotes', 'grRules'})
             Tinfo = strcat(Tinfo,cNameAdd,Cinfo) ;
         % And some are simply seperated by a '|'.
         else
           % Break information from C into parts.
           pipePos = [0 strfind(Cinfo,'|') length(Cinfo)+1] ;
           for iP = 1:length(pipePos)-1
               nowInfo = Cinfo(pipePos(iP)+1:pipePos(iP+1)-1) ; 
               % Only add new pieces of information. 
               if isempty(strfind(Tinfo,nowInfo))  
                 Tinfo = strcat(Tinfo,'|',nowInfo) ;
               end
           end
         end
       % Else, create a new entry. Again note inclusion/absense of
       % model designation. 
       else
           if ismember(rxnFieldNames{iField}, {'subSystems', 'rxnECNumbers', 'rxnReferences', 'rxnNotes', 'grRules'})
               Tinfo = strcat(cNameNew,Cinfo) ;
           else
               Tinfo = Cinfo ;
           end
       end
       % Now update the information in Tmodel with whatever has been made.
       Tmodel.(rxnFieldNames{iField}){tRxn} = Tinfo ; 
     end
   end
   
   % Add information for seperate rxn fields (rev, ub, lb, c)
   for iField = 1:length(sepFieldNames)
       Tmodel.(sepFieldNames{iField}).(cName)(tRxn) = ...
           Cmodel.(sepFieldNames{iField})(cRxn) ; 
   end
   
   % Reverse bounds when reactions have opposite reactants and products.
   % Criteria: reaction has match, reaction is not reversible, and the 
   % products of one match matches the reactants of the other and vice versa.
   if rxnList(cRxn) < tRxns % && ~Cmodel.rev(cRxn)
       % Collect metabolite information from metList assignments.
       cProd = setdiff(metList(Cmodel.S(:,cRxn) > 0) ,ignoreMets) ; 
       tReac = setdiff(find(Tmodel.S(:,tRxn) < 0) , ignoreMets) ;
       cReac = setdiff(metList(Cmodel.S(:,cRxn) < 0) , ignoreMets) ; 
       tProd = setdiff(find(Tmodel.S(:,tRxn) > 0) , ignoreMets);
       if mean(ismember(cProd,tReac)) == 1 && mean(ismember(cReac,tProd)) == 1
           % Do the deed. 
           [Tmodel.lb.(cName)(tRxn),Tmodel.ub.(cName)(tRxn)] = ...
                deal(-Tmodel.ub.(cName)(tRxn), -Tmodel.lb.(cName)(tRxn)) ;
       end
   end
end


%% Add metabolites.
for cMet = 1:cMets
      
    tMet = metList(cMet) ; % Index of met in Tmodel.
    % Update identity array
    Tmodel.Models.(cName).mets(tMet,1) = true ;
    
    % For brand new metabolites.
    if isempty(Tmodel.mets{tMet})
        Tmodel.mets{tMet} = Cmodel.mets{cMet} ;
        Tmodel.metCharge(tMet) = Cmodel.metCharge(cMet) ;
        Tmodel.b(tMet) = Cmodel.b(cMet) ;
        
        % Add in ID info
        Tmodel.metID{tMet} = strcat(Tmodel.metID{tMet},cNameNew, ...
                                    Cmodel.metID{cMet}) ;
    else
        % Add in ID info
        Tmodel.metID{tMet} = strcat(Tmodel.metID{tMet},cNameAdd, ...
                                    Cmodel.metID{cMet}) ;
    end
    
    for iField = 1:length(metFieldNames)
        if iscell(Cmodel.(metFieldNames{iField}))
            % Declare the current info now, makes code cleaner.
            Cinfo = Cmodel.(metFieldNames{iField}){cMet} ;
            Tinfo = Tmodel.(metFieldNames{iField}){tMet} ;
            % Ensure information exists and has not been added
            if ~isempty(Cinfo) && isempty(strfind(Tinfo,Cinfo))
                % Split information up into parts, only adding new info.
                if ~isempty(Tinfo)
                    pipePos = [0 strfind(Cinfo,'|') length(Cinfo)+1] ;
                    for iP = 1:length(pipePos)-1
                        nowInfo = Cinfo(pipePos(iP)+1:pipePos(iP+1)-1) ;
                        if isempty(strfind(Tinfo,nowInfo))
                            Tinfo = strcat(Tinfo,'|',nowInfo) ;
                        end
                    end
                else
                    Tinfo = Cinfo ;
                end
            end     
        % Update information, for each field.
        Tmodel.(metFieldNames{iField}){tMet} = Tinfo ;
        end
    end
    
end

%% ensure must-have fields are there
% Get standard field names. 
fields = TmodelFields ; 
% Reaction related cell field names. 
for ir = 1:length(fields{1})
    if ~isfield(Tmodel,fields{1}{ir})
        Tmodel.(fields{1}{ir}) = cell(size(Tmodel.rxns)) ;
    end
end
% Reaction related double Fields (which are kept model specific).
for ir = 1:length(fields{2})
    if ~isfield(Tmodel,fields{2}{ir})
        Tmodel.(fields{2}{ir}) = zeros(size(Tmodel.rxns)) ;
    end
end
% Metabolite related cell field names. 
for ir = 1:length(fields{3})
    if ~isfield(Tmodel,fields{3}{ir})
        Tmodel.(fields{3}{ir}) = cell(size(Tmodel.mets)) ;
    end
end
% Metabolite related double field names. 
for ir = 1:length(fields{4})
    if ~isfield(Tmodel,fields{4}{ir})
        Tmodel.(fields{4}{ir}) = cell(size(Tmodel.mets)) ;
    end
end

%% output Tmodel
if strcmpi(varargin{1},'NoClean')
    TmodelC = Tmodel ;
else
    % Clean and output Tmodel.
    TmodelC = cleanTmodel(Tmodel) ;
end
