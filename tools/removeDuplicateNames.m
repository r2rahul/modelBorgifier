function infoList = removeDuplicateNames(infoList)
%removeDuplicateNames accepts a cell array of strings, with info seperated
% by a '|' within each cell, and removes the duplicate names in each cell. 
% Also orders the names by size from smallest to largest. Best used after 
% fixNames.
%
%INPUTS
% infoList  Cell array of strings. 
%
%OUPUTS
% infoList
%
%CALLED BY
% cleanTmodel
% verifyModel

%% Get rid of those duplicates!
for iInfo = 1:length(infoList)
    info = infoList{iInfo} ; 
    newInfo = '' ;
    % Break information into parts.
    pipePos = [0 strfind(info,'|') length(info)+1] ;
    newInfoCell = {[]} ; 
    newInfoLength = 0 ; 
    for iP = 1:length(pipePos)-1
        nowInfo = info(pipePos(iP)+1:pipePos(iP+1)-1) ;
        % Check to see if the information exists.
        if isempty(find(strcmp(nowInfo,newInfoCell),1))
            % If not, record the data and it's length
            if isempty(newInfoCell{1})
                newInfoCell{1} = nowInfo ;
                newInfoLength(1) = length(nowInfo) ;
            else
                newInfoCell{length(newInfoCell)+1,1} = nowInfo ;
                newInfoLength(length(newInfoCell),1) = length(nowInfo) ;
            end
        end
    end
    % Sort the data based on length.
    [null lengthOrder] = sort(newInfoLength) ; 
    newInfoCell = newInfoCell(lengthOrder) ; 
    % Assemble cell array into single string. 
    newInfo = newInfoCell{1} ; 
    for iInfoPart = 2:length(newInfoCell) 
        newInfo = [newInfo '|' newInfoCell{iInfoPart}] ;  
    end
    % Replace the information the updated string. 
    infoList{iInfo} = newInfo ; 
end
