% This script finds 
%
%
%INPUTS
%
%
%OUTPUTS
%
%CALLED BY
%




%% Delare variables
load('/modeling/models/Tmodel.mat') 
load('/modeling/models/saccharomyces/imm904_xls+xml.mat')
Cmodel = model ; 

cRxns = length(Cmodel.rxns) ; 
tRxns = length(Tmodel.rxns) ; 

cS = full(Cmodel.S) ; 
tS = full(Tmodel.S) ; 

%% Create matrix of reactions with equal number of involved metabolites.
% Find stoichiometry for reactions in C.
cRxnStoich = cell(cRxns,1) ; 
for iRxn = 1:cRxns
    % Pull out the stoichiometries
    nowStoich = cS(logical(cS(:,iRxn)),iRxn) ;
    % Order each entry for comparison puposes and save in cell.
    cRxnStoich{iRxn} = sort(nowStoich) ; 
end
% Find stoiciometry for reactions in T.
tRxnStoich = cell(tRxns,1) ; 
for iRxn = 1:tRxns
    nowStoich = tS(logical(tS(:,iRxn)),iRxn) ;
    tRxnStoich{iRxn} = sort(nowStoich) ; 
end
clear nowStoich

% Fill in rxnMatch, allotting a 1 where the reactions have the same stoich. 
rxnMatch = zeros(cRxns,tRxns) ; 
for cRxn = 1:cRxns
    for tRxn = 1:tRxns
        % If the two stoics have the same length and their differnce is 0.
        if length(cRxnStoich{cRxn}) == length(tRxnStoich{tRxn})
            if ~(cRxnStoich{cRxn} - tRxnStoich{tRxn})
                % Assign a 1 to that spot.
                rxnMatch(cRxn,tRxn) = 1 ;
            end
        end
    end
end

clear iRxn
%% Compare metabolites
metScores = zeros(size(cS,1),size(tS,1)) ; 
% For each metabolite in C, compare against every metabolite in T. 
for cMet = 1:length(Cmodel.mets)
    for tMet = 1:length(Tmodel.mets)
        % Pull out submatrix of the involved reactions for the 2 mets. 
        subSetRxnMatch = rxnMatch(logical(cS(cMet,:)), ...
                                  logical(tS(tMet,:))) ;
        % Assign a score (= number of reactions with the same
        % stoich/potential number of matches / potenital hits) 
        % Potential hits (number of involved reactions from C or T,
        % whichever is smaller)
        if size(subSetRxnMatch,1) < size(subSetRxnMatch,2)
            matches = sum(logical(sum(subSetRxnMatch,2))) ;
            potentialHits = size(subSetRxnMatch,1) ;
        else
            matches = sum(logical(sum(subSetRxnMatch,1))) ;
            potentialHits = size(subSetRxnMatch,2) ;
        end
        metScores(cMet,tMet) = matches / potentialHits ; 
    end
end

