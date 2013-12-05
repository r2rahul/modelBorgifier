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
function Model = buildTmodel(Model)
%buildTmodel initiates the template model (Tmodel). It adds model
% designations before some of the cell arrays add creates indenty arrays
% that indicate which reactions and metabolties are contained in which
% models. This function should be used after verifyModel. 

%% Declare variables.
modelName = Model.description ;
nameFields = {'rxnID' 'metID' 'subSystems' 'rxnReferences' ...
               'rxnECNumbers' 'grRules'} ;
    
%% Put the model name in for appropriate fields. 
for iF = 1:length(nameFields)
    % Determine which entries are non empty.
    haveData = find(~cellfun(@isempty,Model.(nameFields{iF}))) ;
    Model.(nameFields{iF})(haveData) = ...
        strcat(modelName,':',Model.(nameFields{iF})(haveData)) ; 
end

%% Create indentity structure. 
Model.Models.(modelName).rxns = true(length(Model.rxns),1) ;
Model.Models.(modelName).mets = true(length(Model.mets),1) ;

%% Add genes.
Model.Models.(modelName).genes = Model.genes ; 

%% Move bound info into structures.
rev = Model.rev ;
lb = Model.lb ;
ub = Model.ub ;
c = Model.c ;

Model.rev.(modelName) = rev ;
Model.lb.(modelName) = lb ;
Model.ub.(modelName) = ub ;
Model.c.(modelName) = c ;

%% Make S sparse.
Model.S = sparse(Model.S) ;

%% Remove extra fields that could get in the way.
unwantedFields = {'description' 'genes' 'rxnGeneMat' 'rules'} ; 
Model = rmfield(Model,unwantedFields) ; 