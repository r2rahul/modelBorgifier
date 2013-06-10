function formulaList = fixChemFormulas(formulaList)
%fixChemFormulas formats chemical formulas such that they do not include a
% '1' after elements that only occur once in the molecule. Only operates on
% C, H, N, O, P, S
%
%INPUTS
% formulaList   Cell array of chemical formulas.
%
%OUTPUTS
% formulaList   Same cell array, formated.
%
%CALLED BY      
% verifyModel

%% Find 1's and replace with nothing.
searchString = '(?<=(C|H|N|O|P|S))1(?!\d)' ;
formulaList = regexprep(formulaList,searchString,'') ;




