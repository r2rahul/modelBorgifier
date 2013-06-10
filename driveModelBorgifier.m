% driveAddModel walks through the process of comparing and merging models.
% It is not meant to be used as a function, rather as a guide.



%% Load and verify Cmodel. (The Compare Model). 
fileName = ...
% ['/home/jts/techdata/Bioinformatics/modeling/models' ...
%     '/coli/iAF1260/Ec_iAF1260_flux1.xml'] ;
% fileName = 'T:\Bioinformatics\modeling\models\Mycoplasma pneumoniae\iJW145.xml' ;

% With regular COBRA function.
Cmodel = readCbModel(fileName) ;

% Or with custom written script. 
% Cmodel = readModel_xxxx(fileName) ;

% Verify model is ready for subsequent scripts. 
Cmodel = verifyModel(Cmodel) ;

% Now is a good time to see if this model carries flux. OPTIONAL.


%% Load Tmodel. (The Template Model). 
% Load a matlab workspace with a previously used Tmodel.
if isunix
    load('/models/Tmodel.mat')
elseif ispc
    load('t:\Bioinformatics\modeling\models\Tmodel.mat')
end

% Or use any model as the template model. 
% iBSU = readmodel_ibsu1103 ; 

% If Tmodel is just another model, verify it as well and convert it to a
% proper format for comparison. Also make sure it carries flux. 
if ~isfield(Tmodel,'Models')
    Tmodel = verifyModel(Tmodel) ;
    Tmodel = buildTmodel(Tmodel) ; 
end

%% Compare models. 
% Score Cmodel against Tmodel. 
[Cmodel,Tmodel,score,Stats] = compareCbModels(Cmodel,Tmodel) ;

%% Match models. 
% Initial comparison and matching.
[rxnList,metList,Stats] = reactionCompare(Cmodel,Tmodel,score) ;

% Subsequent comparisons and matching. 
[rxnList,metList,Stats] = reactionCompare(Cmodel,Tmodel,score, ...
                                          rxnList,metList,Stats) ;

%% Merge models and test results.
[TmodelC,Cspawn,Stats1] = evaluateTmodel(Cmodel,Tmodel, ...
                                        rxnList,metList,Stats) ;

%% Extract a model. 
modelToExtract = 'iJW145' ;%'iAF1260' ;
Cspawn = readCbTmodel(modelToExtract,TmodelC) ; 

%% Write to SBML. 
% You will be prompted with a dialoge box for the file name, or you can
% enter it after 'sbml'. This version of writeCbModel does write all the 
% additional information to the .xml.
fileName = 'T:\Bioinformatics\modeling\models\Mycoplasma pneumoniae\iJW145_new.xml' ;
% fileName = ...
% ['/home/jts/techdata/Bioinformatics/modeling/common/cobra' ...
%     '/addModel/saves/iAF1260testoutput.xml'] ;
writeCbModel(Cspawn,'sbml',fileName)

%% Test to see if model can be read back from SBML.
% Note that the readCbModel function doesn't deal with the extra fields 
% that are in the .xml
fileName = ...
['/home/jts/techdata/Bioinformatics/modeling/common/cobra' ...
    '/addModel/saves/iAF1260testoutput.xml'] ;
test = readCbModel(fileName) ;


return
%% save new Tmodel
Tmodel = TmodelC ;
if isunix
    save('/models/Tmodel.mat','Tmodel')
    save(['/models/Tmodel_' datestr(now,'yyyy.mm.dd') '.mat'],'Tmodel')
elseif ispc
    save('t:\Bioinformatics\modeling\models\Tmodel.mat','Tmodel')
    save(['t:\Bioinformatics\modeling\models\Tmodel_' datestr(now,'yyyy.mm.dd') '.mat'],'Tmodel')
end
