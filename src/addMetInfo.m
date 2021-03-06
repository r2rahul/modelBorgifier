function varargout = addMetInfo(varargin)
% ADDMETINFO MATLAB code for addMetInfo.fig
%      ADDMETINFO, by itself, creates a new ADDMETINFO or raises the existing
%      singleton*.
%
%      H = ADDMETINFO returns the handle to a new ADDMETINFO or the handle to
%      the existing singleton*.
%
%      ADDMETINFO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADDMETINFO.M with the given input arguments.
%
%      ADDMETINFO('Property','Value',...) creates a new ADDMETINFO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before addMetInfo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to addMetInfo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help addMetInfo

% Last Modified by GUIDE v2.5 06-Dec-2013 14:19:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @addMetInfo_OpeningFcn, ...
                   'gui_OutputFcn',  @addMetInfo_OutputFcn, ...
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


% --- Executes just before addMetInfo is made visible.
function addMetInfo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to addMetInfo (see VARARGIN)

% Choose default command line output for addMetInfo
handles.output = hObject;

global CMODEL
set(handles.popupmenu_metab,'String',CMODEL.mets)
cfields = fieldnames(CMODEL) ;
for ifn = 1:length(cfields)
    goodfields(ifn) = strncmpi(cfields{ifn},'met',3) ;
end
set(handles.popupmenu_propertylist,'String', cfields(goodfields))

updateeditfield(handles)
set(handles.popupmenu_metab,'Value',varargin{1})
updateeditfield(handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes addMetInfo wait for user response (see UIRESUME)
uiwait(handles.figure1);


function varargout = addMetInfo_OutputFcn(hObject, eventdata, handles) 
% varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_propertylist.
function popupmenu_propertylist_Callback(hObject, eventdata, handles)
updateeditfield(handles)

function popupmenu_propertylist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_prop_Callback(hObject, eventdata, handles)

function edit_prop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton_done_Callback(hObject, eventdata, handles)
close(handles.figure1)

function popupmenu_metab_Callback(hObject, eventdata, handles)
updateeditfield(handles)

function popupmenu_metab_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function updateeditfield(handles)
global CMODEL
nowmet = get(handles.popupmenu_metab,'String') ;
nowmet = nowmet{get(handles.popupmenu_metab,'Value')} ;
nowfield = get(handles.popupmenu_propertylist,'String') ;
nowfield = nowfield{get(handles.popupmenu_propertylist,'Value')} ;

nowval = CMODEL.(nowfield)(strcmp(CMODEL.mets,nowmet)) ;
if isempty(nowval)
    nowval = {''} ;
elseif iscell(nowval)
    nowval = nowval{1} ;
elseif isnumeric(nowval)
    nowval = num2str(nowval) ;
end

set(handles.edit_prop,'String',nowval)



function pushbutton_change_Callback(hObject, eventdata, handles)
global CMODEL
nowmet = get(handles.popupmenu_metab,'String') ;
nowmet = nowmet{get(handles.popupmenu_metab,'Value')} ;
nowfield = get(handles.popupmenu_propertylist,'String') ;
nowfield = nowfield{get(handles.popupmenu_propertylist,'Value')} ;

nowval = get(handles.edit_prop,'String') ;

if iscell(CMODEL.(nowfield)(strcmp(CMODEL.mets,nowmet)))
    nowval = {nowval} ;
elseif isnumeric(CMODEL.(nowfield)(strcmp(CMODEL.mets,nowmet)))
    nowval = str2double(nowval) ;
end

CMODEL.(nowfield)(strcmp(CMODEL.mets,nowmet)) = nowval;

