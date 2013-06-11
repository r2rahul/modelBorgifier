function cNum = countC(nowFormula)
%countC finds number of carbons in a formula, if any.
%
% cNum = countC(nowFormula)

% Remove alternate formulas.
nowFormula = regexprep(nowFormula,'/|,*','') ;

% Look for C's
cNum = regexp(nowFormula,'C\d*','match') ;
if isempty(cNum)
    cNum = 0 ;
else
    if length(cNum{1}) == 1 ;
        cNum = 1 ;
    else
        cNum = str2double(cNum{1}(2:end)) ;
    end
end


% cpos = strfind(nowFormula,'C') ;
% if isempty(cpos)
%     cNum = 0 ;
% else
%     for ic = 1:length(cpos)
%         if sum(strcmp(nowFormula(cpos(ic)+1), {'a' 'o' 'r' 'u' 'd' 's' 'l'})) == 0 % if not Ca, Co, Cr, Cu, Cd, Cs, Cr Cl
%             if isempty(str2num(nowFormula(cpos(ic)+1))) % 1 C
%                 cNum(ic) = 1 ;
%             elseif isempty(str2num(nowFormula(cpos(ic)+1:cpos(ic)+2))) % 2 to 9 C
%                 cNum(ic) = str2num(nowFormula(cpos(ic)+1)) ;
%             elseif isempty(str2num(nowFormula(cpos(ic)+1:cpos(ic)+3))) % 10 to 99 C
%                 cNum(ic) = str2num(nowFormula(cpos(ic)+1:cpos(ic)+2)) ;
%             else % 100 to 999 C
%                 cNum(ic) = str2num(nowFormula(cpos(ic)+1:cpos(ic)+3)) ;
%             end
%         else
%             cNum(ic) = 0 ;
%         end
%     end
% end
% cNum = sum(cNum) ;
% clear ic ; clear nowformula ; clear cpos ;