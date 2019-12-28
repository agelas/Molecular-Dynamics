function varargout = mdGUI(varargin)
% MDGUI MATLAB code for mdGUI.fig
%      MDGUI, by itself, creates a new MDGUI or raises the existing
%      singleton*.
%
%      H = MDGUI returns the handle to a new MDGUI or the handle to
%      the existing singleton*.
%
%      MDGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MDGUI.M with the given input arguments.
%
%      MDGUI('Property','Value',...) creates a new MDGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mdGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mdGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mdGUI

% Last Modified by GUIDE v2.5 03-May-2019 16:03:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mdGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @mdGUI_OutputFcn, ...
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


% --- Executes just before mdGUI is made visible.
function mdGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mdGUI (see VARARGIN)

% Choose default command line output for mdGUI
handles.output = hObject;

% Update handles structure
handles.go=0;
guidata(hObject, handles);

% UIWAIT makes mdGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mdGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    handles.res=md(36,6,0.2);
    handles.res.draw;
    %axes(handles.axes1)
    ax=handles.axes1;
