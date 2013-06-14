% driveModelBorgifier walks through the process of comparing and merging 
% models. It is not meant to be used as a function, rather as a guide.
% Please reference the manual and the help information in the individual 
% scripts for additional information.

%% Load and verify Cmodel. (The Compare Model).
% Load with regular COBRA function.
if isunix
    fileName = '/test/iBSU1103.xml';
else ispc
    fileName = '\test\iBSU1103.xml';
end
Cmodel = readCbModel(fileName); 

% Or with custom written script. 
Cmodel = readModel_xxxx(fileName);

% Verify model is ready for subsequent scripts. 
Cmodel = verifyModel(Cmodel);

% If model has SEED IDs, use the databases to fill in information.
rxnFileName = '/addModelgit/test/SEED_db/ModelSEED-reactions-db.csv'; 
cpdFileName = '/addModelgit/test/SEED_db/ModelSEED-compounds-db.csv';
Cmodel = addSEEDInfo(Cmodel,rxnFileName,cpdFileName); 
 
% Now is a good time to see if this model carries flux. OPTIONAL.
solverOkay = changeCobraSolver('glpk','LP'); 
CmodelSol = optimizeCbModel(Cmodel); 

%% Load Tmodel. (The Template Model). 
% Load a matlab workspace with a previously used Tmodel.
if isunix
    load('/Tmodel.mat')
elseif ispc
    load('\Tmodel.mat')
end

% Or use any model as the template model.
if isunix
    fileName = '/test/iJO1366.xml';
else ispc
    fileName = '\test\iJO1366.xml';
end

% Load with regular COBRA function.
Tmodel = readCbModel(fileName);

% If Tmodel is just another model, verify it as well and convert it to a
% proper format for comparison. Also make sure it carries flux. 
if ~isfield(Tmodel,'Models')
    Tmodel = verifyModel(Tmodel);
    TmodelSol = optimizeCbModel(Tmodel); 
    Tmodel = buildTmodel(Tmodel); 
end

%% Compare models. 
% Score Cmodel against Tmodel. This can taken a few hours. 
[Cmodel,Tmodel,score,Stats] = compareCbModels(Cmodel,Tmodel);

%% Match models.
% Initial comparison and matching.
[rxnList,metList,Stats] = reactionCompare(Cmodel,Tmodel,score);

% OPTIONAL. Declare mets from Cmodel with comps not in Tmodel as new.
metList1 = newCompsNewMets(metList,Cmodel,Tmodel);

% Subsequent comparisons and matching. 
[rxnList,metList,Stats] = reactionCompare(Cmodel,Tmodel,score, ...
                                          rxnList,metList,Stats);

%% Merge models and test results.
[TmodelC,Cspawn,Stats1] = evaluateTmodel(Cmodel,Tmodel, ...
                                         rxnList,metList,Stats);

%% Extract a model. 
modelToExtract = 'iJO1366';
Cspawn = readCbTmodel(modelToExtract,TmodelC); 

%% Write to SBML. 
% You will be prompted with a dialoge box for the file name, or you can
% enter it after 'sbml'. This version of writeCbModel does write all the 
% additional information to the .xml.

fileName = 'model.xml';
writeCbModel(Cspawn,'sbml')

%% Test to see if model can be read back from SBML.
% Note that the readCbModel function doesn't deal with the extra fields 
% that are in the .xml
fileName = '/saves/testoutput.xml'];
test = readCbModel(fileName);

%% save new Tmodel
Tmodel = TmodelC ;
if isunix
    save('/models/Tmodel.mat','Tmodel')
    save(['/models/Tmodel_' datestr(now,'yyyy.mm.dd') '.mat'],'Tmodel')
elseif ispc
    save('t:\Bioinformatics\modeling\models\Tmodel.mat','Tmodel')
    save(['t:\Bioinformatics\modeling\models\Tmodel_' datestr(now,'yyyy.mm.dd') '.mat'],'Tmodel')
end
