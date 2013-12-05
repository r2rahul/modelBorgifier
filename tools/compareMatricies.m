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
function FluxCompare = compareMatricies(CMODEL,Cspawn)
% Creates a structure with the matricies before and after merging.
%
% Called by mergeModels
%

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

% find protons and water
CMODELprotonpos = logical(strncmpi(FluxCompare.CmetsSort,'h[',2) + ...
    strncmpi(FluxCompare.CmetsSort,'h+[',3) + ...
    strncmpi(FluxCompare.CmetsSort,'c0065[',6) + ...
    strncmpi(FluxCompare.CmetsSort,'proton[',7) ) ;

CMODELwaterpos = logical(strncmpi(FluxCompare.CmetsSort,'h2o[',4) + ...
    strncmpi(FluxCompare.CmetsSort,'water[',6) ) ;

ignorepos = logical(CMODELprotonpos + CMODELwaterpos) ;
% ignore protons and water
FluxCompare.diffS(ignorepos, :) = 0 ;

end