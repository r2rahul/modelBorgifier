function fields = TmodelFields(~)
% TmodelFields returns the field names of the Tmodel structure. This is
% to help allievate control between different scripts.
%
% [rxnFields,numField,metFields,allFields] = TmodelFields(~)
%
%INPUTS
% 
%OUTPUTS
% Fields        A column cell array that contains the following in order. 
% rxnFields
% rNumField
% metFields
% mNumFields
% allFields
%
%CALLED BY
% verifyModel

% Declare fields
fields = cell(4,1) ;

% Reaction related fields that are cell arrays.
rxnFields = {'rxns' 'rxnID' 'rxnNames' 'subSystems' 'rxnECNumbers' ...
             'rxnKEGGID' 'rxnSEEDID' 'rxnEquations' ...
             'rxnReferences' 'rxnNotes' 'grRules'}' ;
fields{1,1} = rxnFields ; 
             
% Reaction related fields that are numeric arrays. 
rNumFields = {'rev' 'lb' 'ub' 'c'}' ;
fields{2,1} = rNumFields ; 
             
% Metabolite related fields that are cell arrays
metFields = {'mets' 'metID' 'metNames' 'metFormulas' ...
             'metKEGGID' 'metSEEDID' ...
             'metChEBIID' 'metPubChemID' 'metInChIString'}' ;
fields{3,1} = metFields ;

% Metabolite related fields that are numeric arrays.
mNumFields = {'metCharge'} ; 
fields{4,1} = mNumFields ; 

% All field names in correct order (27 fields).
allFields = {'rxns' 'mets' 'S' 'rev' 'lb' 'ub' 'c' ...
             'rxnID' 'rxnNames' 'subSystems' 'rxnEquations' ...
             'rxnECNumbers' 'rxnKEGGID' 'rxnSEEDID' ...
             'rxnReferences' 'rxnNotes' 'grRules' ...
             'metID' 'metNames' 'metFormulas' 'metCharge' ...
             'metKEGGID' 'metSEEDID' 'metChEBIID' 'metPubChemID' ...
             'metInChIString' ...
             'genes' 'description'}' ;
fields{5,1} = allFields ; 
                
