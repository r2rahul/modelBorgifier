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
function Model = orderModelFields(Model)
% orderModelFields puts the fields in the correct order. 
%
%INPUTS
% Model
%
%OUTPUTS
% Model
%
%CALLS
% TmodelFields
%
%CALLED BY      
% verifyModel

%% Order fields.
% Get correct order. 
fields = TmodelFields ; 
correctOrder = fields{5} ;

% Current Field names in the model. 
fieldNames = fieldnames(Model) ;

% Extra fields should go at the end.
fieldOrder = zeros(length(fieldNames),1) ;
extraFieldIndex = length(correctOrder) + 1 ;

% Determine current field position
for iField = 1:length(fieldNames)
    fieldIndex = strcmp(fieldNames{iField},correctOrder) ;
    fieldIndex = find(fieldIndex) ;
    if fieldIndex
        fieldOrder(fieldIndex) = iField ;
    else
        fieldOrder(extraFieldIndex) = iField ;
        extraFieldIndex = extraFieldIndex + 1 ;
    end
end

% fix vFieldOrderS
fieldOrder(fieldOrder == 0) = [] ;

% Organize
Model = orderfields(Model,fieldOrder) ;
