function [Model] = addSEEDInfo(Model) 
% This function takes the compound and reaction database files from SEED
% and compares them to a model, adding any additional information it can
% find and checking to see if current information is in agreement with the
% databases. It returns an augmented model. 
%
% 
%INPUTS
% Tmodel
%
%OUTPUTS
% Tmodel
% Stats     Structure with relevant stats on what has been updated. 






%% Declare variables.
% Compound database filename. 
cpdFileName = ('/modeling/models/seed/ModelSEED-compound-db.csv') ; 
% Reaction database filename. 
rxnFileName = ('/modeling/models/seed/ModelSEED-reaction-db.csv') ; 





%% Compare Tmodel metabolites to metabolite database.

% Check metabolites in T with a SEED ID, and then see if these metabolites
% already have the information available from the database. Add info if it
% does not, throw flag if the information is different. 


% Now use KEGGIDs and an anchor point for moving in information.
% Exclude metabolites that have SEEDIDs (these should have been convered
% above)


% Repeat, using names (not sure which comparison, or series of comparisons,
% is most appropriate here).


% Replace compounds in Tmodel that use the SEED cpdXXXXX name with the 
% abbreviations from the database. 



%% Compare Tmodel reactions to reaction database. 

% Check reactions from Tmodel that have SEEDIDs, fill in information from
% database if it is not already present. Flag what is differnt. 


% Do the same, but using KEGGIDs as the anchor. 


% Replace the reactions in Tmodel which use the SEED rxnXXXXX name with the
% abbreviatios from the database. 

