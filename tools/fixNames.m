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
function nameList = fixNames(nameList)
%fixNames standardizes names, ie, met, rxn and long names for them as well.
% fixNames makes the names lowercase, removes obtuse characters and
% whitespace and replaces them with underscores, and then removes starting
% and trailing underscores.
%
%INPUTS
% nameList      Cell array of names.
%
%OUTPUTS
% nameList      Same cell array but with fixed names. 
%
%CALLED BY
% verifyModel

%% Fix those names!
% Make everything lowercase. 
nameList = lower(nameList) ;

% Remove wierd characters.
nameList = regexprep(nameList,' |-(?!($|\||\[))|,|:|\(|\)','_') ;

% Remove [ and ] that aren't part of the compartment lable.
% ie ] that are not at the end of the word.
nameList = regexprep(nameList,'\](?!($|\|))','_') ;
% And [ that don't look ahead do a letter and a ].
nameList = regexprep(nameList,'\[(?!.\]$)','_') ;

% Consolodate resulting underscores.
nameList = regexprep(nameList,'_____|____|___|__','_') ;

% Delete ones at beginning and end of names. 
nameList = regexprep(nameList, '^_|\|_|_\||_$|_(?=\[\w\]$)','') ;

