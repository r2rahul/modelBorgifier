function [TmodelC,Cspawn,Stats] = evaluateTmodel(CmodelIn,TmodelIn, ...
                                                 rxnList,metList,Stats)
% evaluateTmodel checks Tmodel for duplicate reactions and other mistakes,
% that may have occured during reaction and metabolite matching. It
% resolves these problems and merges the models, and confirms that Cmodel
% is the same after being removed from the merged model. It also provides 
% some statistics on the merging process and resulting combined model.
%
% [TmodelC,Cspawn,Stats] = evaluateTmodel(CmodelIn,TmodelIn, ...
%                                         rxnList,metList,Stats)
%
%INPUT
% Cmodel    
% Tmodel
% rxnList
% metList
% Stats     Structure that comes from reactionCompare. Weighting
%           information can be used and additional information addended. 
%
%OUTPUT
% TmodelC   Combined C and Tmodel. 
% Cspawn    Cmodel after it has been removed from the TmodelC
% Stats     Structure of information regarding the merging and also stats
%           on the combined model.
%
%CALLS
% TmodelStats

%% Declare variables. 
global CMODEL TMODEL
CMODEL = CmodelIn ;
TMODEL = TmodelIn ; 

% Allow pausing for some of the UI. 
pause on

%% Check that metList is okay. 
% Make sure all mets have been declared and that 2 mets from C aren't
% matched to the same met in T. Continue until everything is a okay.
reviewMets = 1 ; 
while reviewMets 
    % Find mets that have not been declared.
    reviewMets = find(metList == 0) ;

    % Find mets in C that have been matched to the same metabolite in T. 
    uniqMetIndex = unique(metList) ; 
    count = histc(metList,uniqMetIndex) ;
    dups = uniqMetIndex(count>1) ;
    for iDup = 1:length(dups)
        metIndex = find(metList == dups(iDup)) ;  
        reviewMets(end+1:length(reviewMets)+length(metIndex)) = metIndex ; 
    end

    % Review all above mets with GUI. 
    if ~isempty(reviewMets)
        fprintf('Problems within metList, resolve with GUI.\n')
        fprintf('Press the any key to continue.\n')
        pause
        RxnInfo.rxnIndex = 0 ; % Tells GUI just look at mets.  
        RxnInfo.rxnList = rxnList ;
        RxnInfo.metList = metList ;
        RxnInfo.metIndex = reviewMets ; 
        metList = metCompare(RxnInfo) ;
    end
end

% Force good numbering for new metabolites
TmetNum = length(TMODEL.mets) ;
metList(metList > TmetNum) = TmetNum+(1:sum(metList>TmetNum)) ;

%% Reconsider 


%% Reactions that don't have the same stoich between Cmodel and Cspawn
checkSimilarity = 1 ; 
while checkSimilarity
    % Combine models and create spawn of Cmodel. 
    [TmodelC,Cspawn] = mergeAndSpawn(CMODEL,TMODEL,rxnList,metList) ;
    
    FluxCompare = compareMatricies ;
    
    % If there are differences pause and let the user know what is up. 
    if sum(sum(FluxCompare.diffS))
        % Plots the differences between the S matricies. 
        figure
        subplot(2,2,1)
        spy(FluxCompare.CmodelS)
        title([CMODEL.description ' before merging.'])
        subplot(2,2,2)
        spy(FluxCompare.CspawnS)
        title([CMODEL.description ' after merging.']) 
        subplot(2,1,2)
        spy(FluxCompare.diffS)
        title('Difference.') 

        % Pause and tell the users what is happening.
        fprintf('Difference between the matricies. Reactions from C\n')
        fprintf('that do not have the same stoich before and after\n')
        fprintf('merging will be declared as new and readded.\n')
        fprintf('Press the any key to continue.\n')
        pause

        % Find the wrong reactions, mark them as new in RxnList, and loop.
        diffStoichs = logical(sum(abs(FluxCompare.diffS),1)) ;
        rxnList(FluxCompare.CrxnsSorti(diffStoichs)) = 0 ; 
    else
        % Set the flag no not check for differences again.
        fprintf('Matricies are now equal before and after merging.\n')
        % Add rxn and metList to Stats.
        Stats.rxnList = rxnList ; 
        Stats.metList = metList ; 
        % Turn off flag. 
        checkSimilarity = 0 ;
    end
end
    
    % Creates a structure with the matricies before and after merging. 
    function FluxCompare = compareMatricies() 
    clear so
    % Sort rxns and mets so they are in the same order.
    [FluxCompare.CrxnsSort,FluxCompare.CrxnsSorti] = sort(CMODEL.rxnID) ; 
    [FluxCompare.SrxnsSort,FluxCompare.SrxnsSorti] = sort(Cspawn.rxnID) ;
    [FluxCompare.CmetsSort,FluxCompare.CmetsSorti] = sort(CMODEL.metID) ; 
    [FluxCompare.SmetsSort,FluxCompare.SmetsSorti] = sort(Cspawn.metID) ;

    % Create sorted S matricies for Cmodel before and after it was added. 
    FluxCompare.CmodelS = CMODEL.S(FluxCompare.CmetsSorti, ...
                                   FluxCompare.CrxnsSorti) ;
    FluxCompare.CspawnS = Cspawn.S(FluxCompare.SmetsSorti, ...
                                   FluxCompare.SrxnsSorti) ;

    % Find the differences between the matricies. 
    FluxCompare.diffS = abs(FluxCompare.CmodelS) - ...
                        abs(FluxCompare.CspawnS) ;
    end

%% Checking flux calculations.

% Are objective functions the same?
if ~(find(CMODEL.c(FluxCompare.CrxnsSorti)) == ...
        find(Cspawn.c(FluxCompare.SrxnsSorti)))
    fprintf('Objective function index not conserved after merging\n')
end

% Are bounds the same?
lbDiff = find(CMODEL.lb(FluxCompare.CrxnsSorti) - ...
              Cspawn.lb(FluxCompare.SrxnsSorti)) ;
ubDiff = find(CMODEL.ub(FluxCompare.CrxnsSorti) - ...
              Cspawn.ub(FluxCompare.SrxnsSorti)) ;

solverOK = changeCobraSolver('glpk','LP') ; 
if ~solverOK
    fprintf('LP solver is not set correctly\n')
end

% Testing flux that comes out.
FBACmodel = optimizeCbModel(CMODEL,'max') ; 
FBACspawn = optimizeCbModel(Cspawn,'max') ;

% Order the flux values so they can be easily compared. 
FBACmodel.xSort = FBACmodel.x(FluxCompare.CrxnsSorti) ; 
FBACspawn.xSort = FBACspawn.x(FluxCompare.SrxnsSorti) ;

% Difference in fluxes. 
FluxCompare.fluxDiff = abs(FBACmodel.xSort - FBACspawn.xSort) ;
FluxCompare.fluxDiff(FluxCompare.fluxDiff < 1e-7) = 0  ;

% If the flux values 
if (FBACmodel.f ~= FBACspawn.f) || ~isempty(find(FluxCompare.fluxDiff,1))
    fprintf('Fluxes not the same. Check Stats.FBA\n')
end

Stats.FluxCompare = FluxCompare ; 
Stats.FBACmodel = FBACmodel ; 
Stats.FBACspawn = FBACspawn ; 

%% Organize and Get stats on combined Tmodel.
TmodelC = organizeModelCool(TmodelC) ;
Stats = TmodelStats(TmodelC,Stats) ;
end

%% Subfunctions.
% Combines C and Tmodel and pulls C back out, returns the merge and spawn.
function [TmodelC,Cspawn] = mergeAndSpawn(CMODEL,TMODEL,rxnList,metList)
    % Combine the models. 
    TmodelC = addToTmodel(CMODEL,TMODEL,rxnList,metList) ;

    % Create Cmodel again from the combined model. 
    Cspawn = readCbTmodel(CMODEL.description,TmodelC) ;  
end




% %% Checking for duplicate reactions.
% % Check for reactions with the exact same stoichiometry.
% identicalStoich = ones(nRxns) ; 
% for iRxn = 1:nRxns-1
%     for jRxn = iRxn+1:nRxns
%         identicalStoich(iRxn,jRxn) = sum(abs(TmodelC.S(:,iRxn) - ...
%                                              TmodelC.S(:,jRxn))) ;  
%     end
% end
% identicalStoich = sparse(~logical(identicalStoich)) ;
% 
% figure
% spy(identicalStoich1)
% title('Reactions with identical stoichiometry.')
% 
% [Stats.sameStoichX,Stats.sameStoichY] = find(identicalStoich) ;
%     
% 
% % Compare based on SEED ID. 
% 
% 
% % Check based on reweighted data. 