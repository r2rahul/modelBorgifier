function varargout = BSSsvm(varargin)
% BSSSVM MATLAB code for BSSsvm.fig
%      BSSSVM, by itself, creates a new BSSSVM or raises the existing
%      singleton*.
%
%      H = BSSSVM returns the handle to a new BSSSVM or the handle to
%      the existing singleton*.
%
%      BSSSVM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BSSSVM.M with the given input arguments.
%
%      BSSSVM('Property','Value',...) creates a new BSSSVM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BSSsvm_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BSSsvm_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BSSsvm

% Last Modified by GUIDE v2.5 28-Feb-2012 09:56:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BSSsvm_OpeningFcn, ...
                   'gui_OutputFcn',  @BSSsvm_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before BSSsvm is made visible.
function BSSsvm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BSSsvm (see VARARGIN)

% Choose default command line output for BSSsvm
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
%%  logo
%load the background image into Matlab
%if image is not in the same directory as the GUI files, you must use the
%full path name of the iamge file
logoImage = imread('logo_holmes120.png');
%select the axes
axes(handles.axes_logo);
%place image onto the axes
image(logoImage);
%remove the axis tick marks
axis off

global d
global fm
% get list of fieldnames from fm
filterlist = fieldnames(fm) ;
% check if fields contain a single number
goodlist = [] ;
for il = 1:length(filterlist) ;
    if isnumeric(fm(1).(filterlist{il})) && length(fm(1).(filterlist{il})) == 1
        goodlist(end+1) = il ;
    end
end
% write fields containing a single number to popup menu
filterlist = filterlist(goodlist) ;
set(handles.listbox_featurelist,'String',filterlist) ;

clear filterlist
% set list of selection
filterlist1 = fieldnames(d.selection) ;
filterlist2 = fieldnames(fm(1).selection) ;
for il = 1:length(filterlist1)
    filterlist{il,1} = ['mrHolmes_' filterlist1{il}] ;
end
for il2 = 1:length(filterlist2)
    filterlist{il+il2,1} = ['drWatson_' filterlist2{il2}] ;
end

set(handles.popupmenu_grouplist1,'String',filterlist) ;
set(handles.popupmenu_grouplist2,'String',['Alle anderen' ; filterlist] ) ;

set(handles.edit_svmname,'String',['svm_' datestr(now,'yymmdd')])
uiwait

    
% UIWAIT makes BSSsvm wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = BSSsvm_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout = {[]} ;

function edit_paramqc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function listbox_featurelist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_grouplist1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_grouplist2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_gamma_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_cost_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_svmname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_svmname_Callback(hObject, eventdata, handles)
function edit_cost_Callback(hObject, eventdata, handles)
function edit_gamma_Callback(hObject, eventdata, handles)
function edit_paramqc_Callback(hObject, eventdata, handles)
function popupmenu_grouplist1_Callback(hObject, eventdata, handles)
function popupmenu_grouplist2_Callback(hObject, eventdata, handles)
function listbox_featurelist_Callback(hObject, eventdata, handles)


function pushbutton_cancel_Callback(hObject, eventdata, handles)
close(handles.figure1)

function pushbutton_go_Callback(hObject, eventdata, handles)
try
    % get d, lf, and fm from global workspace
    global d
    global fm
    global lf
    set(handles.text_busy,'Visible','on')
    pause(0.0001)
    
    [data, groups, grouplabels] = svmprepare(d, fm, handles) ;
    
    % train model on single measurements (fm) because the dispersity of
    % replicates is needed for the training
    model = svmtrain(grouplabels, data(:,groups)',  [' -c ' num2str(2^str2num(get(handles.edit_cost,'String'))) ' -g ' num2str(2^str2num(get(handles.edit_gamma,'String')))]) ;
    
    % make prediction
    [prediction, accuracy, decision_value] = svmpredict(ones(size(data,2),1), data', model) ;
    nowfieldname = get(handles.edit_svmname,'String') ;
    for ic = 1:d.n
        nowmean = mean(decision_value(d.pos{ic})) ;
        d.([nowfieldname '_dv'])(ic,1) = nowmean ;
        if nowmean > 0
            d.selection.(nowfieldname)(ic,1) = 1 ;
        else
            d.selection.(nowfieldname)(ic,1) = 0 ;
        end
    end
    
    featurelist = get(handles.listbox_featurelist,'String') ;
    featurelist = {featurelist{get(handles.listbox_featurelist,'Value')'}}' ;
    grouplist1 = get(handles.popupmenu_grouplist1,'String') ;
    grouplist1 = grouplist1{get(handles.popupmenu_grouplist1,'Value')'} ;
    grouplist2 = get(handles.popupmenu_grouplist2,'String') ;
    grouplist2 = grouplist2{get(handles.popupmenu_grouplist2,'Value')'} ;

    flout = [] ;
    for ifl = 1:length(featurelist)
        flout = [flout featurelist{ifl} '; '] ;
    end
    flout = flout(1:end-2) ;
    
    lf{end+1,1} = [datestr(now,'yyyy.mm.dd HH:MM') ': SVM berechnet mit Features; ' flout ] ;
    lf{end+1,1} = ['Gruppe 1: ' grouplist1 ', Gruppe2: ' grouplist2] ;
    
    set(handles.text_busy,'Visible','off')
    close(handles.figure1)
catch exception
    BSSerrormessage(exception)
end

function pushbutton_testparam_Callback(hObject, eventdata, handles)
try
    % get d from global workspace
    set(handles.text_busy,'Visible','on')
    pause(0.0001)
    global d
    global fm
    
    [data, groups, grouplabels] = svmprepare(d, fm, handles) ;
    
    qc = svmtrain(grouplabels, data(:,groups)',  ['-v 5 -c ' num2str(2^str2num(get(handles.edit_cost,'String'))) ' -g ' num2str(2^str2num(get(handles.edit_gamma,'String')))]) ;
    
    set(handles.edit_paramqc,'String', num2str(qc, '%g'))
    set(handles.text_busy,'Visible','off')
catch exception
    BSSerrormessage(exception)
end

function pushbutton_autoparameters_Callback(hObject, eventdata, handles)
try
    % get d from global workspace
    global d
    global fm
    set(handles.text_busy,'Visible','on')
    pause(0.000001)
    
    [data, groups, grouplabels] = svmprepare(d, fm, handles) ;
    
    stepsize = 10 ;
    costrange = [-30 : stepsize: 30] ;
    gammarange = [-20 : stepsize: 40] ;
    lastqc = 0 ;
    bestqc = 0 ;
    improve = true ;
    backupcnt = 0 ;
    while improve
        improve = false ;
        backupcnt = backupcnt + 1 ;
        clear nowqc ; clear qc ;
        ic = 0 ;
        for cost = costrange
            ic = ic + 1 ;
            ig = 0 ;
            for gamma = gammarange
                ig = ig + 1 ;
                nowqc(ic,ig) = svmtrain(grouplabels, data(:,groups)',  ['-v 3 -c ' num2str(2^cost) ' -g ' num2str(2^gamma) ' -q']) ;
                if nowqc(ic,ig) > bestqc
                    bestcost = cost ;
                    bestgamma = gamma ;
                    set(handles.edit_paramqc,'String',num2str(bestqc)) ;
                end
            end
        end
        
        for ic2 = 2:ic-1
            for ig2 = 2:ig-1
                qc(ic2-1 , ig2-1) = mean(mean(nowqc([ic2-1 : ic2+1],[ig2-1:ig2+1]))) ;
            end
        end
        bestqc2 = max(max(qc)) ;
        bestpos = qc == bestqc2 ;
        bestpos = [round(mean(find(sum(bestpos,2)))) round(mean(find(sum(bestpos,1)))) ] + 1 ;
        
        %     figure
        %     surf(gammarange, costrange, nowqc)
        
        stepsize = stepsize * 0.5
        costrange = [costrange(bestpos(1)-1) : stepsize : costrange(bestpos(1)+1)]
        gammarange = [gammarange(bestpos(2)-1) : stepsize : gammarange(bestpos(2)+1)]
        
        if lastqc < bestqc2 - 0.03
            improve = true ;
            lastqc = bestqc ;
        end
        if backupcnt > 5
            improve = false ;
        end
    end
    
    set(handles.edit_cost,'String',num2str(bestcost)) ;
    set(handles.edit_gamma,'String', num2str(bestgamma)) ;
    
    set(handles.text_busy,'Visible','off')
catch exception
    BSSerrormessage(exception)
end


function [data, groups, grouplabels] = svmprepare(d,fm, handles)
% get selected features
featurelist = get(handles.listbox_featurelist,'String') ;
featurelist = {featurelist{get(handles.listbox_featurelist,'Value')'}}' ;

% compile feature matrix
for ifl = 1:length(featurelist)
    data(ifl,:) = [fm(:).(featurelist{ifl})] ;
end ; clear ifl ;

% normalize 
for id = 1:size(data,1)
    % set minimum to zero
    data(id,:) = data(id,:) - min(data(id,:)) ;
    % if data is far from normally distributed then use log of data rather
    % than data
    if abs(log2(  mean(data(id,:)) / median(data(id,:))  )) > 1
        data(id,:) = data(id,:) + 1 ;
        data(id,:) = log10(data(id,:)) ;
    end
    % set maximum to one
    data(id,:) = data(id,:) ./ max(data(id,:)) ;
end ; clear id ;

% get selected group for selection 1
grouplist = get(handles.popupmenu_grouplist1,'String') ;
grouplist = grouplist{get(handles.popupmenu_grouplist1,'Value')'} ;

if strfind(grouplist,'mrHolmes') == 1
    groups = [d.pos{d.selection.(grouplist(10:end)) ~= 0} ]'  ;
    grouplabels = ones(length([d.pos{d.selection.(grouplist(10:end)) ~= 0} ]  ),1) ;
elseif strfind(grouplist,'drWatson') == 1
    for iw = 1:length(fm)
        refpos(iw) = (fm(iw).selection.(grouplist(10:end)) ~= 0) ;
    end
    groups = find(refpos)' ;
    grouplabels = ones(length(find(refpos)),1) ;
end

% get selected group for selection 2
grouplist = get(handles.popupmenu_grouplist2,'String') ;
grouplist = grouplist{get(handles.popupmenu_grouplist2,'Value')'} ;

if strfind(grouplist,'mrHolmes') == 1
    groups = [groups ; [d.pos{d.selection.(grouplist(10:end)) ~= 0} ]' ] ;
    grouplabels = [grouplabels ; 2*ones(length([d.pos{d.selection.(grouplist(10:end)) ~= 0} ]  ),1) ];
    
elseif strfind(grouplist,'drWatson') == 1
    for iw = 1:length(fm)
        refpos(iw) = (fm(iw).selection.(grouplist(10:end)) ~= 0) ;
    end
    groups = [groups ; find(refpos)' ];
    grouplabels = [grouplabels ; 2*ones(length(find(refpos)),1) ] ;
    
elseif strcmp(grouplist,'Alle anderen')
    if length(groups) > (length(fm) / 3)
        grouplabels = [grouplabels ; zeros(length(fm) - length(groups) ,1)] ;
        groups = [groups ; setdiff( (1:length(fm)), groups)'  ] ;
    else
        grouplabels = [grouplabels ; zeros(2*length(groups) ,1)] ;
        groups = [groups ; randsample(setdiff( (1:length(fm)), groups), 2*length(groups))'  ] ;
    end
end
