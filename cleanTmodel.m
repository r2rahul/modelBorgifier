function Tmodel = cleanTmodel(Tmodel) 
%cleanTmodel Reorganizes and checks Tmodel for completeness. To be used 
% immediately after addToTmodel. 
%
%Tmodel = cleanTmodel(Tmodel) 
%
%INPUTS
% Tmodel
%
%OUTPUTS
% Tmodel
%
%CALLED BY
% addToTmodel
%
%CALLS
% TmodelFields
% buildRxnEquations
% removeDuplicateNames
% makeNamesUnique

%% Declare variables.
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

% Model names.
moNames = fieldnames(Tmodel.Models) ;

% Number of reactions and number of metabolites.
nRxns = length(Tmodel.rxns) ;
nMets = length(Tmodel.mets) ;


%% Lengthen rev, lb, ub, and c arrays for old models.
for iField = 1:length(rNumFields)
    for iMo = 1:length(moNames)
        nowArray = Tmodel.(rNumFields{iField}).(moNames{iMo}) ;
        if length(nowArray) ~= nRxns
            nowArray(end:nRxns) = 0 ;
            Tmodel.(rNumFields{iField}).(moNames{iMo}) = nowArray;
        end
    end
end

%% Lengthen reaction and metabolite identity arrays, convert to logicals.
for iMo = 1:length(moNames)
    % Reactions.
    rxnsArray = Tmodel.Models.(moNames{iMo}).rxns ;
    if length(rxnsArray) ~= nRxns
        rxnsArray(end+1:nRxns) = 0 ;
    end
    Tmodel.Models.(moNames{iMo}).rxns = logical(rxnsArray) ;
    
    % Metabolites. 
    metsArray = Tmodel.Models.(moNames{iMo}).mets ;
    if length(metsArray) ~= nMets
        metsArray(end+1:nMets) = 0 ;
    end
    Tmodel.Models.(moNames{iMo}).mets = logical(metsArray) ;
end

%% Make sure all names are unique before reorganizing.
% Rebuild reaction equations so there is a hint for renaming rxns. 
Tmodel = buildRxnEquations(Tmodel) ; 
fprintf('Checking if reaction IDs (.rxns) are unique.\n') ; 
Tmodel.rxns = makeNamesUnique(Tmodel.rxns,Tmodel.rxnEquations) ;
fprintf('Checking if metabolite IDs (.mets) are unique.\n') ;
Tmodel.mets = makeNamesUnique(Tmodel.mets,Tmodel.metNames) ; 

% Remove duplicate met IDs.
Tmodel.metID = removeDuplicateNames(Tmodel.metID) ; 

% Rebuild reaction equations. 
Tmodel = buildRxnEquations(Tmodel) ;

%% Reorder reactions and metabolites alphabetically.
% Reorder reaction related lists.
[~,rxnI] = sort(Tmodel.rxns) ;
for iField = 1:length(rxnFields)
    Tmodel.(rxnFields{iField}) = Tmodel.(rxnFields{iField})(rxnI) ;
end
for iField = 1:length(rNumFields)
    for iMo = 1:length(moNames)
        Tmodel.(rNumFields{iField}).(moNames{iMo}) = ...
            Tmodel.(rNumFields{iField}).(moNames{iMo})(rxnI) ;
    end
end

% Reorder metabolite related lists.
[~,metI] = sort(Tmodel.mets) ;
for iField = 1:length(metFields)
    Tmodel.(metFields{iField}) = Tmodel.(metFields{iField})(metI) ;
end
for iField = 1:length(mNumFields)
    Tmodel.(mNumFields{iField}) = ...
        Tmodel.(mNumFields{iField})(metI) ;
end

% Reorder identity arrays.
for iMo = 1:length(moNames)
    Tmodel.Models.(moNames{iMo}).rxns = ... 
        Tmodel.Models.(moNames{iMo}).rxns(rxnI) ;
    Tmodel.Models.(moNames{iMo}).mets = ...
        Tmodel.Models.(moNames{iMo}).mets(metI) ;
end

% Reorder S matrix.
Tmodel.S = Tmodel.S(:,rxnI) ;
Tmodel.S = Tmodel.S(metI,:) ;

%% Share information between the same metabolite in different compartments.
% Metabolite fields to share information between
share = {'metNames' 'metFormulas' 'metKEGGID' 'metSEEDID' 'metChEBIID' ...
         'metPubChemID' 'metInChIString'} ;
          
metsNoComp = cell(nMets,1) ;
for iMet = 1:nMets
    metsNoComp{iMet} = Tmodel.mets{iMet}(1:end-3) ;
end
uniqMets = unique(metsNoComp,'first') ;

% Create combined formation fields for each metabolite and redistribute.
for iMet = 1:length(uniqMets)
  % Find sister metabolites
  sMets = find(strcmp(uniqMets{iMet},metsNoComp)) ;
  if length(sMets) > 1
    % Create combined information for each field.
    for iF = 1:length(share)
      info = {''} ;
      for iSis = 1:length(sMets)
        % If the field is not empty, break information into parts
        % (as seperated by pipes), and add that information in if
        % it is not already present. 
        if ~isempty(Tmodel.(share{iF}){sMets(iSis)})
          pipePos = [0 ...
                     strfind(Tmodel.(share{iF}){sMets(iSis)},'|') ...
                     length(Tmodel.(share{iF}){sMets(iSis)})+1 ] ;
          for j = 1:length(pipePos)-1
            nowInfo = Tmodel.(share{iF}){sMets(iSis)} ...
                      (pipePos(j)+1:pipePos(j+1)-1) ;
            % If the info does not already exist.
            if isempty(strfind(info{1},nowInfo))
              % If it is the first piece of information or additional. 
              if isempty(info{1})
                info{1} = nowInfo ;
              else
                info{1} = strcat(info{1}, '|', nowInfo) ;
              end
            end
          end
        end
      end
      for iSis = 1:length(sMets)
          Tmodel.(share{iF}){sMets(iSis)} = info{1} ;
      end
      clear info
    end
  end
end

% Now remove duplicate information which may have snuck in. 
for iField = 1:length(share)
    Tmodel.(share{iField}) = removeDuplicateNames(Tmodel.(share{iField})) ; 
end
  
%% Remove duplicate reaction longnames and IDs
Tmodel.rxnNames = removeDuplicateNames(Tmodel.rxnNames) ;
Tmodel.rxnID = removeDuplicateNames(Tmodel.rxnID) ;

%% Clear leftover fields from the matching process
deleteFields = {'rxnComp' 'metNums' 'rxnMetNames'} ; 
existFields = isfield(Tmodel,deleteFields) ;
for iF = 1:length(existFields)
    if existFields(iF)
        Tmodel = rmfield(Tmodel,deleteFields{iF}) ;
    end
end


%% Unused code (but helpful for cleaning up models if there are mistakes.
% It is a little messy, sorry

% %%
% 
% % Fixing empty rxnID
% noRxnID = cellfun('isempty',Tmodel.rxnID) ;
% noRxnID = find(noRxnID) ;
% noRxnIDName = Tmodel.rxnNames(noRxnID) ;
% for i = 1:length(noRxnIDName)
%     noRxnIDName{i} = noRxnIDName{i}(10:end) ;
% end
% clear i
% for i = 1:length(noRxnID)
%     test(i,1) = strcat('iBSU1103:',Tmodel.rxns(noRxnID(i))) ;
% end
% Tmodel.rxnID(noRxnID) = test;
% noRxnIDrxns = Tmodel.rxns(noRxnID) ;
% 
% 
% 
% metPosInT = zeros(length(noRxnID),1) ;
% for i = 1:length(noRxnID)
%     metPosInT(i) = find(Tmodel.S(:,noRxnID(i))) ;
% end
% metNamesInT = Tmodel.mets(metPosInT) ;
% test = unique(metNamesInT) ;
% 
% 
% % Bad mets (ones that do not have handedness)
% for i = 1:length(noRxnIDName)
%     rxnIndexBSU(i) = find(strcmp(noRxnIDName{i},iBSU1103.rxnNames)) ;
% end
% rxnNameBSU = rxnNameBSU(:) ;
% badMets = find(strcmp(noRxnIDName{i},iBSU1103.rxnNames)) ;
% badMetEx = find(strcmp('ala[e] Export',noRxnIDName)) ;
% iBSU1103.rxnNames(badMets)
% metNamesInT(badMetEx)
% Tmodel.rxnNames(noRxnID(badMetEx))
% test = unique(noRxnIDName) ;
% test2 = iBSU1103.rxnNames(rxnNameBSU) ;
% test3 = iBSU1103.rxnEquations(rxnNameBSU) ;
% iBSU1103.rxnNames(1640)
% 
% find(strcmp('glu[e] Export',noRxnIDName))
% gluExRxns = find(strcmp(noRxnIDName{i},iBSU1103.rxnNames)) ;
% gluExEq = iBSU1103.rxnEquations(gluExRxns) ;
% for j = 1:2
%     gluExMetPos(j) = find(iBSU1103.S(:,gluExRxns(j)))
% end
% gluExMetNames = iBSU1103.metNames{1085} ;
% 
% 
% % No rxn name
% noRxnName = cellfun('isempty',Tmodel.rxnNames) ;
% noRxnName = find(noRxnName) ;
% noRxnNameID = Tmodel.rxnID(noRxnName) ;
% noRxnNameEq = Tmodel.rxnEquations(noRxnName) ;
% for i = 1:length(noRxnName)
%     test(i,1) = strcat('iBSU1103:',Tmodel.rxns(noRxnName(i))) ;
% end
% Tmodel.rxnNames(noRxnName) = test ;
% 
% for i = 1:length(test)
%     rxnNameBSU(i) = find(strcmp(test{i},iBSU1103.rxnID)) ;
% end
% rxnNameBSU = rxnNameBSU(:) ;
% test2 = iBSU1103.rxnNames(rxnNameBSU) ;
% test3 = iBSU1103.rxnEquations(rxnNameBSU) 
% 
% test = Tmodel.rxnNames(noRxnss) ;
% 
% % No subsystems
% noRxnss = cellfun('isempty',Tmodel.subSystems) ;
% noRxnss = find(noRxnss) ;
% 
% % No Formulas.
% noMetForm = cellfun('isempty',Tmodel.metFormulas) ;
% noMetForm = find(noMetForm) ; 
% test = Tmodel.mets(noMetForm) ;
% 
% 
% %% Unused Code
% 
% % sharedRxnsLog(Tmodel.Models.iAF1260.rxns == Tmodel.Models.iBSU1103.rxns) = 1;
% % sharedRxnsLog = logical(sharedRxnsLog(:)) ;
% % sharedRxnIndexs = find(sharedRxnsLog) ;
% % sharedRxns = Tmodel.rxns(sharedRxnsLog) ;
% % sharedRxnsNo = sum(sharedRxnsLog) ;
% % uniqRxnsLog = logical(uniqRxnsLog(:)) ;
% % uniqRxnIndexs = find(sharedRxnsLog) ;
% % uniqRxns = Tmodel.rxns(uniqRxnsLog) ;
% % uniqRxnsNo = sum(uniqRxnsLog) ;
% 
% % Find reactions not in either model.
% % uhoh = zeros(nRxns,1) ;
% % for i = 1:length(sharedRxns)
% %     if Tmodel.Models.iAF1260.rxns(i) == 0 && ...
% %             Tmodel.Models.iBSU1103.rxns(i) == 0 
% %          uhoh(i) = 1 ;
% %     end
% % end
% % uhoh = find(uhoh) ;
% 
% % Find mets not in either model. 
% uhoh = zeros(nMets,1) ;
% for i = 1:length(shared)
%     if Tmodel.Models.iAF1260.mets(i) == 0 && ...
%             Tmodel.Models.iBSU1103.mets(i) == 0 
%          uhoh(i) = 1 ;
%     end
% end
% uhoh = find(uhoh) ;
% uhoh(11) = 2307 ;
% uhoh = uhoh(1:10)
% 
% Tmodel.mets(uhoh)
% Tmodel.metID(uhoh)
% Tmodel.Models.iBSU1103.mets(uhoh) = 1 ;