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