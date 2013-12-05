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
function [TmodelC,Cspawn,Stats, CMODEL] = mergeModels(CmodelIn,TmodelIn, ...
    rxnList,metList,Stats,score)
% mergeModels checks Tmodel for duplicate reactions and other mistakes,
% that may have occured during reaction and metabolite matching. It
% resolves these problems and merges the models, and confirms that Cmodel
% is the same after being removed from the merged model. It also provides
% some statistics on the merging process and resulting combined model.
%
% [TmodelC,Cspawn,Stats] = mergeModels(CmodelIn,TmodelIn, ...
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
% addToTmodel
% readCbTmodel
% TmodelStats
% optimizeCbModel

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
CMODELoriginal = CMODEL ;
try
checkSimilarity = 1 ;
while checkSimilarity
    % Combine models and create spawn of Cmodel.
    [TmodelC,Cspawn,CMODEL] = mergeAndSpawn(CMODEL,TMODEL,rxnList,metList) ;
    
    FluxCompare = compareMatricies(CMODELoriginal,Cspawn) ;
    
    % compare bounds
    boundDiff = [] ;
    for ir = 1:length(CMODEL.rxns)
        relativeDirection = FluxCompare.CmodelS(:,ir) ./ FluxCompare.CspawnS(:,ir) ;
        nowmets = setdiff(find(FluxCompare.CmodelS(:,ir)) ,...
            find(strncmp( Cspawn.mets(FluxCompare.SmetsSorti),'h[',2)) ); % ignore protons
        relativeDirection = relativeDirection(nowmets)  ;
        if ~isempty(nowmets)
            if mean(relativeDirection) == 1
                if ~CMODEL.lb(FluxCompare.CrxnsSorti(ir)) == Cspawn.lb(FluxCompare.SrxnsSorti(ir)) || ...
                        ~CMODEL.ub(FluxCompare.CrxnsSorti(ir)) == Cspawn.ub(FluxCompare.SrxnsSorti(ir))
                    boundDiff(end+1) = ir ;
                end
            elseif mean(relativeDirection) == -1
                if ~CMODEL.lb(FluxCompare.CrxnsSorti(ir)) == -Cspawn.ub(FluxCompare.SrxnsSorti(ir)) || ...
                        ~CMODEL.ub(FluxCompare.CrxnsSorti(ir)) == -Cspawn.lb(FluxCompare.SrxnsSorti(ir))
                    boundDiff(end+1) = ir ;
                end
            else
                rxnList(ir) = -1 ;
                metList(nowmets(abs(relativeDirection)~=1)) = 0 ;
            end
        end
    end
    if ~isempty(boundDiff)
        fprintf('bounds not conserved after merging\nErrors in: \n' )
        rxnList(boundDiff) = -1 ;
    end
    
    % If there are differences pause and let the user know what is up.
    if sum(sum(FluxCompare.diffS)) || ~isempty(boundDiff) || (sum(rxnList == -1) > 0)
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
        fprintf('merging need to be reviewed again.\n')
        fprintf('Press the any key to continue.\n')
        pause
        
        % Find the wrong reactions and metabolies,
        % mark them as undeclared, call GUI, and loop.
        
        diffmets = logical(sum(abs(FluxCompare.diffS),2)) ;
        metList(FluxCompare.CmetsSorti(diffmets)) = 0 ;
        diffrxns= logical(sum(abs(FluxCompare.diffS),1) + sum(FluxCompare.CmodelS(diffmets,:)) + sum(FluxCompare.CspawnS(diffmets,:)) ) ;
        rxnList(FluxCompare.CrxnsSorti(diffrxns)) = -1 ;
        [rxnList, metList, Stats] = reactionCompare(CMODEL,TMODEL,score, rxnList, metList, Stats) ;
        
    else
        % Set the flag no not check for differences again.
        fprintf('Matricies are now equal before and after merging.\n')
        % Turn off flag.
        checkSimilarity = 0 ;
    end
    % Add rxn and metList to Stats.
    Stats.rxnList = rxnList ;
    Stats.metList = metList ;
end
catch % return savely if re-matching is aborted by user
    return
end

%% Checking flux calculations.
% Are objective functions the same?
if ~(find(CMODEL.c(FluxCompare.CrxnsSorti)) == ...
        find(Cspawn.c(FluxCompare.SrxnsSorti)))
    fprintf('Objective function index not conserved after merging\n')
end

% % % % Are bounds the same?
% % % problemspersist = true ;
% % % while problemspersist
% % %     [TmodelC,Cspawn,CMODEL] = mergeAndSpawn(CMODEL,TMODEL,rxnList,metList) ;
% % %     problemspersist = false ;
% % %     boundDiff = [] ;
% % %     for ir = 1:length(CMODEL.rxns)
% % %         relativeDirection = CMODEL.S(FluxCompare.CmetsSorti,FluxCompare.CrxnsSorti(ir)) ./ ...
% % %             Cspawn.S(FluxCompare.SmetsSorti,FluxCompare.SrxnsSorti(ir)) ;
% % %         nowmets = setdiff( find(CMODEL.S(FluxCompare.CmetsSorti,FluxCompare.CrxnsSorti(ir))) ,...
% % %             find(strncmp( Cspawn.mets,'h[',2)) ); % ignore protons
% % %         relativeDirection = relativeDirection(nowmets)  ;
% % %         if ~isempty(nowmets)
% % %             if mean(relativeDirection) == 1
% % %                 if ~CMODEL.lb(FluxCompare.CrxnsSorti(ir)) == Cspawn.lb(FluxCompare.SrxnsSorti(ir)) || ...
% % %                         ~CMODEL.ub(FluxCompare.CrxnsSorti(ir)) == Cspawn.ub(FluxCompare.SrxnsSorti(ir))
% % %                     boundDiff(end+1) = ir ;
% % %                 end
% % %             elseif mean(relativeDirection) == -1
% % %                 if ~CMODEL.lb(FluxCompare.CrxnsSorti(ir)) == -Cspawn.ub(FluxCompare.SrxnsSorti(ir)) || ...
% % %                         ~CMODEL.ub(FluxCompare.CrxnsSorti(ir)) == -Cspawn.lb(FluxCompare.SrxnsSorti(ir))
% % %                     boundDiff(end+1) = ir ;
% % %                 end
% % %             else
% % %                 disp('This message should not appear. It means that there might be an error in the Stoichiometry after all...')
% % %                 rxnList(ir) = -1 ;
% % %                 metList(nowmets(abs(relativeDirection)~=1)) = 0 ;
% % %                 problemspersist = true ;
% % %             end
% % %         end
% % %     end
% % %     if ~isempty(boundDiff)
% % %         fprintf('bounds not conserved after merging\nErrors in: \n' )
% % %         rxnList(boundDiff) = 0 ;
% % %         disp(FluxCompare.CrxnsSort(boundDiff))
% % %         problemspersist = true ;
% % %     end
% % %     [rxnList, metList, Stats] = reactionCompare(CMODEL,TMODEL,score, rxnList, metList, Stats) ;
% % %     % Add rxn and metList to Stats.
% % %     Stats.rxnList = rxnList ;
% % %     Stats.metList = metList ;
% % % end

solverOK = changeCobraSolver('glpk','LP') ;
if ~solverOK
    fprintf('LP solver is not set correctly\n')
end

% Testing flux that comes out.
try
    FBACmodel = optimizeCbModel(CMODEL,'max','one') ;
    FBACspawn = optimizeCbModel(Cspawn,'max','one') ;

    % Order the flux values so they can be easily compared.
    FBACmodel.xSort = FBACmodel.x(FluxCompare.CrxnsSorti) ;
    FBACspawn.xSort = FBACspawn.x(FluxCompare.SrxnsSorti) ;

    % Difference in fluxes.
    FluxCompare.fluxDiff = abs(FBACmodel.xSort - FBACspawn.xSort) ;
    FluxCompare.fluxDiff(FluxCompare.fluxDiff < 1e-7) = 0  ;
catch
    return
end

% If the flux values
if (FBACmodel.f ~= FBACspawn.f) || ~isempty(find(FluxCompare.fluxDiff,1))
    fprintf('Fluxes not the same. Check Stats.FBAoverview\n')
    [~, csort] = sort(abs(FBACmodel.x),'descend') ;
    [~, ssort] = sort(abs(FBACspawn.x),'descend') ;
    Stats.FBAoverview = [CMODEL.rxns(csort) CMODEL.rxnEquations(csort) num2cell(FBACmodel.x(csort)) ...
        Cspawn.rxns(ssort) Cspawn.rxnEquations(ssort) num2cell(FBACspawn.x(ssort)) ] ;
end

Stats.FluxCompare = FluxCompare ;
Stats.FBACmodel = FBACmodel ;
Stats.FBACspawn = FBACspawn ;

%% Organize and Get stats on combined Tmodel.
TmodelC = cleanTmodel(TmodelC) ;
try
    TmodelC = organizeModelCool(TmodelC) ;
catch
    disp('TmodelC too large for cool organization of stoichiometric matrix')
end
Cspawn = readCbTmodel(CMODEL.description,TmodelC) ;
Stats = TmodelStats(TmodelC,Stats) ;
end

%% Subfunctions.
% Combines C and Tmodel and pulls C back out, returns the merge and spawn.
function [TmodelC,Cspawn,CMODEL] = mergeAndSpawn(CMODEL,TMODEL,rxnList,metList)
% Combine the models.
[TmodelC, CMODEL] = addToTmodel(CMODEL,TMODEL,rxnList,metList,'NoClean') ;

% Create Cmodel again from the combined model.
Cspawn = readCbTmodel(CMODEL.description,TmodelC) ;
end
