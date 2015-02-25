function varargout = spm_pc(varargin)

% SPM_PC M-file for spm_pc.fig
% the Patient Classification toolbox is designed
% to classify one or several subjects in comparison
% with a control group - bootstrapped CI are computed
% for the control group and subjects are classified
%
% the toolbox requires SPM5 and the Matlab statistic (optional)
% and image processing (required) toolbox
% cyril pernet 26/06/2008


%% initialize variables
global defaults 
spm_defaults; ;

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spm_pc_OpeningFcn, ...
                   'gui_OutputFcn',  @spm_pc_OutputFcn, ...
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

% --- Executes just before spm_pc is made visible.
function spm_pc_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% define handles used for the save callback
handles.one_sample_ci         = 0;
handles.ci_images             = 0;
handles.path                  = pwd;
handles.estimate_H0           = 0;
handles.classification_images = 0;
handles.controls              = 0;
handles.patients              = 0;
handles.alpha_value           = 0;
handles.Nboot                 = 0;
guidata(hObject, handles);

function varargout = spm_pc_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%% --- Executes on button press in spm_select.
function spm_select_Callback(hObject, eventdata, handles)
global defaults  

handles.controls = spm_select(Inf,'.*\.img$','Select Images to compute CI from (controls)');
handles.patients = spm_select(Inf,'.*\.img$','Select Images to classify (patients)');

if isempty(handles.controls) == 1
    handles.controls = 0;
end

if isempty(handles.patients) == 1
    handles.patients = 0;
end
    
guidata(hObject, handles);


%% --- Executes on button press in investigate.
function investigate_Callback(hObject, eventdata, handles)
global defaults  
spmpc_ci_test(handles.path, handles.controls)
guidata(hObject, handles);


%% --- Executes on button press in bootstrap_ci.
function bootstrap_ci_Callback(hObject, eventdata, handles)
global defaults  
if handles.controls == 0
    handles.controls = spm_select(Inf,'.*\.img$','Select Images to compute CI from (controls)');
end
[handles.alpha_value, handles.Nboot] = spmpc_ci(handles.one_sample_ci, handles.ci_images, handles.controls, handles.path, handles.alpha_value, handles.Nboot, 0);
guidata(hObject, handles);


% --- Executes on button press in one_sample_ci.
function one_sample_ci = one_sample_ci_Callback(hObject, eventdata, handles)
global defaults  
one_sample_ci = get(hObject,'Value'); % returns toggle state of one_sample_ci
if one_sample_ci == 1
    disp('one sample ci computation turned on');
    handles.one_sample_ci = 1;
else
    disp('one sample ci computation turned off');
    handles.one_sample_ci = 0;
end
guidata(hObject, handles);

% --- Executes on button press in ci_images.
function ci_images = ci_images_Callback(hObject, eventdata, handles)
global defaults 
ci_images = get(hObject,'Value'); % returns toggle state of one_sample_ci
if ci_images == 1
    disp('generation of ci images turned on');
    handles.ci_images = 1;
else
    disp('generation of ci images turned off');
    handles.ci_images = 0;
end
guidata(hObject, handles);


%%  Executes on button press in classify.
function classify_Callback(hObject, eventdata, handles)
global defaults 
if handles.patients == 0
    handles.patients = spm_select(Inf,'.*\.img$','Select Images to compute CI from');
end
spmpc_classification(handles.estimate_H0, handles.classification_images, handles.controls, handles.patients, handles.path, handles.alpha_value, handles.Nboot);
guidata(hObject, handles);

% --- Executes on button press in estimate_h0.
function estimate_h0_Callback(hObject, eventdata, handles)
global defaults 
estimate_H0 = get(hObject,'Value'); % returns toggle state of one_sample_ci
if estimate_H0 == 1
    disp('generation of ci images turned on');
    disp('beware this takes a very long time');
    handles.estimate_H0 = 1;
else
    disp('generation of ci images turned off');
    handles.estimate_H0 = 0;
end
guidata(hObject, handles);

% --- Executes on button press in classification_images.
function classification_images_Callback(hObject, eventdata, handles)
global defaults 
classification_images = get(hObject,'Value'); % returns toggle state of one_sample_ci
if classification_images == 1
    disp('generation of ci images turned on');
    handles.classification_images = 1;
else
    disp('generation of ci images turned off');
    handles.classification_images = 0;
end
guidata(hObject, handles);

%% --- Executes on button press in cluster_threshold.
function cluster_threshold_Callback(hObject, eventdata, handles)
coordinates = spmpc_cluster(1);
save cluster_coordinates.txt coordinates -ascii
guidata(hObject, handles);


%% --- Executes on button press in cluster_inference.
function cluster_inference_Callback(hObject, eventdata, handles)
max_cluster_size = spmpc_cluster(2);
save max_cluster_size_H0.txt max_cluster_size -ascii
guidata(hObject, handles);


% --- Executes on button press in cross_validation.
function cross_validation_Callback(hObject, eventdata, handles)
spmpc_K_fold(handles.path, handles.controls, handles.patients, handles.alpha_value, handles.Nboot)
guidata(hObject, handles);


%% Executes on button press in cd.
function cd_Callback(hObject, eventdata, handles)
global defaults
handles.path = uigetdir(pwd,'CD');
cd (handles.path);
guidata(hObject, handles);


%% Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
cd (handles.path);
web(['file://' which('spmpc_help.html')]);
guidata(hObject, handles);


%% Executes on button press in quit.
function quit_Callback(hObject, eventdata, handles)
clc; disp('bye ...');
uiresume; delete(handles.figure1)
close all










