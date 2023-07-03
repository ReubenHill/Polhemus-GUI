function varargout = DIGIGUI(varargin)
% -------------------------------------------------------------------------
% Copyright (C) Reuben W. Nixon-Hill (formerly Reuben W. Hill)
%
% -------------------------------------------------------------------------
%                   -- DIGIGUI v1.3.0 for Matlab R2016b  --
% -------------------------------------------------------------------------
%
% For the Polhemus PATRIOT digitiser, attached to stylus pen with button.
%
% A list of points to digitise is imported from a text file, where each
% point is on a new line.
%
% The baud rate is set via the variable "BaudRate" in
% DIGIGUI_OutputFcn and has default value 115200.
%
% Points are digitised by pressing the stylus button.
%
% The default atlas requires 5 reference cardinal points: 'Nasion','Inion',
% 'Ar','Al' and 'Cz'. Once measured, an Atlas reference baby head, with
% these points marked, is mapped onto the graph display of points.
%
% Before the allignment of the head model, a coordinate transform is done:
% 1: place the 'inion' at the origin
% 2: rotate the 'Al' into the y axis
% 3: rotate the 'Ar' into the xy plane about the new 'Inion'-'Al' y axis
% 4: Rotate about the 'Al'-'Ar' axis to bring the 'Nasion' into the xy
%    plane, thus alligning the inion and nasion.
% This coordinate transform is then applied to all measured points.
%
% At present this is the only atlas that has been tested - see the
% 'default_atlas' directory for what is needed for a different atlas. I
% expect that some work will need to be done to make different atlasses
% importable when the program has been compiled: in particular changes to
% the working directory to try and include new .m files.
%
% Once a list of points has been found, the atlas can be realigned at any
% time. This is useful if one wishes to remeasure an atlas point.
%
% A list of points in absolute coordinates of the Polhemus can be specified
% to measure as 'expected coordinates'. These are specified in a CSV file
% as a list of locations and x-y-z coordinates in centimetres - see
% 'expected_coordinates_example.txt' for the format required. The tolerance
% for accepting a measurement is specified by providing a location called
% "Tolerance" with the X coordinate used as the tolerance. The list of
% points must correspond to at least some of the imported locations,
% otherwise the expected coordinates will not be imported and an error will
% be displayed.
%
% A CSV file, mat file, or excel spreadsheet of points and their XYZ coords is
% outputted to a location of the user's choosing.
%
%
% MATLAB GUIDE Generated comments:
% DIGIGUI MATLAB code for DIGIGUI.fig
%      DIGIGUI, by itself, creates a new DIGIGUI or raises the existing
%      singleton*.
%
%      H = DIGIGUI returns the handle to a new DIGIGUI or the handle to
%      the existing singleton*.
%
%      DIGIGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIGIGUI.M with the given input arguments.
%
%      DIGIGUI('Property','Value',...) creates a new DIGIGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DIGIGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DIGIGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DIGIGUI

% Last Modified by GUIDE v2.5 03-Jul-2023 15:48:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DIGIGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @DIGIGUI_OutputFcn, ...
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


% --- Executes just before DIGIGUI is made visible.
function DIGIGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DIGIGUI (see VARARGIN)

% Choose default command line output for DIGIGUI
handles.output = hObject;

disp(['Working dir: ', pwd]);
disp(['ctfroot: ', ctfroot]);

%-------------------Get the executable/.m directory--------------------

if isdeployed
    [status, result] = system('path');
    handles.currentDir = char(regexpi(result, 'Path=(.*?);', ...
        'tokens', 'once'));
end

%-------------------Get the default user location--------------------

if ispc
    handles.userDir = getenv('USERPROFILE');
else
    handles.userDir = getenv('HOME');
end

%--------------------define close request function----------------------
%the function "CloseFcn" that I define now runs when quitting the gui
set(gcf,'CloseRequestFcn',{@CloseFcn,handles});

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DIGIGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = DIGIGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%-------------------DISABLE HEAD ALIGNMENT BUTTON----------------------
set(handles.HeadAlign,'Enable','off');

%-------------------CHECK FOR EXISTING SERIAL OBJECT----------------------

% Look for any existing serial port objects and warn user they may be
% deleted
button = 'OK';
if(~isempty(instrfindall))
    %show message box
    msgline1 = 'The Polhemus Patriot device will now be looked for on available COM ports.';
    msgline2 = 'Warning: existing serial port objects in MATLAB will be deleted if running as a MATLAB script.';
    msgline3 = 'Press OK to continue';
    msg = sprintf('%s\n\n%s\n\n%s',msgline1,msgline2,msgline3);
    disp(msgline2);
    button = questdlg(msg,'Initialising...','OK','Cancel','OK');
end

if(~strcmp(button,'OK'))
    %Quit the program
    disp('Now quitting by user request.');
    guidata(hObject, handles);
    CloseFcn(hObject,eventdata,handles);
    return;
end

%------------------------CREATE SERIAL OBJECT---------------------------

% Create serial object and set baud rate
BaudRate = 115200;

% Find interface objects that are set to 'on' i.e. enabled...
InterfaceObj=findobj(handles.figure1,'Enable','on');
% ... and turn them off.
set(InterfaceObj,'Enable','off');

% find serial com port
[handles.COMport, handles.sensors] = FindPatriotSerial(BaudRate);

% Re-enable the interface objects.
set(InterfaceObj,'Enable','on');

if(handles.COMport ~= 0) %patriot found
    handles.serial = serial(handles.COMport,'BaudRate', BaudRate);
else

    %-------------------QUIT & ERROR IF DEVICE NOT FOUND--------------------
    str1 = 'Polhemus Patriot Device not found or communicated with successfully.';
    str2 = ['Check the device is on and its baud rate is set to '...
        sprintf('%i',BaudRate) ...
        ' on both the hardware switches of the device and the settings of'...
        ' any USB link cable used.'];
    str3 = ['If running in MATLAB, try restarting MATLAB to scan for new'...
        ' serial devices.'];
    str4 = ['Also consider turning the device off and on. Take care to'...
        ' give the device time to reinitialise before trying again.'];
    errstr = sprintf('%s\n\n%s\n\n%s\n\n%s',str1,str2,str3,str4);
    % display error message
    uiwait(errordlg(errstr,'Polhemus Communications Initialisation Error'));
    % quit if com port not found
    guidata(hObject, handles);
    CloseFcn(hObject,eventdata,handles);
    return;

end


%--------------------INITIALISE HANDLES VARIABLES--------------------
% Set the initial point count to 0. This is incremented before each
% measured head point until the last head point is measured.
handles.point_count = 0;
handles.point_count_history = [];

% this is true when opening save dialogues for example
handles.disable_measurements = false;

% this is true if the locations list has been edited.
handles.editedLocationsList = false;
% this is true if the atlas point names has been edited.
handles.editedAtlasPoints = false;

handles.error_distance = 0.1;
handles.double_tap_error_enabled = false;
if ~isdeployed
    atlas_dir = fullfile(pwd, 'default_atlas');
else
    % Note that for some reason the default atlas directory is
    % located directly in ctfroot
    atlas_dir = fullfile(ctfroot, 'default_atlas');
end

%Get settings
try
    if ~isdeployed
        settings_loc = fullfile(pwd, 'settings.mat');
    else
        settings_loc = fullfile(ctfroot, 'DIGIGUI', 'settings.mat');
    end
    disp(['looking for settings.mat at ', settings_loc]);
    load(settings_loc,'error_distance','double_tap_error_enabled','atlas_dir','save_expected_coords', 'only_measure_with_expected_coords', 'insert_by_location_name', 'previousStartNumber', 'previousEndNumber', 'previousStartLetter', 'previousEndLetter');
    disp(['error_distance: ', num2str(error_distance)]);
    disp(['double_tap_error_enabled: ', num2str(double_tap_error_enabled)]);
    disp(['atlas_dir: ', atlas_dir]);
    disp(['save_expected_coords: ', num2str(save_expected_coords)]);
    disp(['only_measure_with_expected_coords: ', num2str(only_measure_with_expected_coords)]);
    disp(['insert_by_location_name: ', insert_by_location_name]);
    disp(['previousStartNumber: ', previousStartNumber]);
    disp(['previousEndNumber: ', previousEndNumber]);
    disp(['previousStartLetter: ', previousStartLetter]);
    disp(['previousEndLetter: ', previousEndLetter]);
catch
    uiwait(errordlg('Settings missing or not found, falling back on default settings.','Settings Error'));
    try
        load('default_settings.mat','error_distance','double_tap_error_enabled','atlas_dir','save_expected_coords','only_measure_with_expected_coords', 'insert_by_location_name', 'previousStartNumber', 'previousEndNumber', 'previousStartLetter', 'previousEndLetter');
        disp(['error_distance: ', num2str(error_distance)]);
        disp(['double_tap_error_enabled: ', num2str(double_tap_error_enabled)]);
        if ~isdeployed
            atlas_dir = fullfile(pwd, atlas_dir);
        else
            atlas_dir = fullfile(ctfroot, atlas_dir);
        end
        disp(['atlas_dir: ', atlas_dir]);
        disp(['save_expected_coords: ', num2str(save_expected_coords)]);
        disp(['only_measure_with_expected_coords: ', num2str(only_measure_with_expected_coords)]);
        disp(['insert_by_location_name: ', insert_by_location_name]);
        disp(['previousStartNumber: ', previousStartNumber]);
        disp(['previousEndNumber: ', previousEndNumber]);
        disp(['previousStartLetter: ', previousStartLetter]);
        disp(['previousEndLetter: ', previousEndLetter]);
    catch
        uiwait(errordlg('Default settings missing or not found! Quitting.','Settings Error'));
        %Quit the gui
        guidata(hObject, handles);
        CloseFcn(hObject,eventdata,handles);
        return
    end
end

handles.error_distance = error_distance;
handles.double_tap_error_enabled = double_tap_error_enabled;
handles.save_expected_coords = save_expected_coords;
handles.only_measure_with_expected_coords = only_measure_with_expected_coords;
handles.menu_options_insert_by_location_name.Checked = insert_by_location_name;
if strcmp(insert_by_location_name, 'on')
    handles.menu_options_insert_by_number.Checked = 'off';
else
    handles.menu_options_insert_by_number.Checked = 'on';
end
handles.previousStartNumber = previousStartNumber;
handles.previousEndNumber = previousEndNumber;
handles.previousStartLetter = previousStartLetter;
handles.previousEndLetter = previousEndLetter;

%Import atlas
noatlas = true;
isdefault = false;
while noatlas
    try
        %load other data needed for headpoint plotting - these are the required
        %names after getting atlas_dir from settings.mat
        [landmarks, landmark_names, mesh] = load_atlas_data(atlas_dir);
        if good_atlas_data(landmarks, landmark_names, mesh)
            handles.atlas_dir = atlas_dir;
            handles.AtlasLandmarks = landmarks;
            handles.AtlasLandmarkNames = landmark_names;
            handles.mesh = mesh;
            noatlas = false;
        end
    catch
        try
            disp('Bad atlas directory, using default');
            load(which('default_settings.mat'), 'atlas_dir');
            if ~isdeployed
                atlas_dir = fullfile(pwd, atlas_dir);
            else
                atlas_dir = fullfile(ctfroot, atlas_dir);
            end
            disp(['atlas_dir: ', atlas_dir]);
            disp(['is default atlas dir? ' num2str(isdefault)]);
            assert(~isdefault);
            isdefault = true;
        catch
            uiwait(errordlg('Default atlas directory missing or not found! Quitting.','Settings Error'));
            %Quit the gui
            guidata(hObject, handles);
            CloseFcn(hObject,eventdata,handles);
            return
        end
    end
end

% Include the atlas directory to get necessary functions
addpath(handles.atlas_dir);

%--------------------HEADPOINTS TO DIGITISE INPUT-----------------------

if ~isdeployed
    saved_location_names_loc = fullfile(pwd, 'savedLocationNames.mat');
else
    saved_location_names_loc = fullfile(ctfroot, 'DIGIGUI', 'savedLocationNames.mat');
end
try
    disp(['looking for savedLocationNames.mat at ', saved_location_names_loc])
    load(saved_location_names_loc, 'locations');
catch
    uiwait(warndlg('Could not find previously used location list.',...
        'Location Warning','modal'));

    % Ask user to load a location list file. Note that the default deployed
    % location is handles.currentDir instead of handles.userDir since the
    % example locations list file can be found in the executable directory
    % (handles.currentDir).
    if ~isdeployed
        [filename,pathname] = ...
            uigetfile({'*.txt;*.dat;*.csv', ...
            'Text Files (*.txt) (*.dat) (*.csv)'} ...
            ,['Select Location List File - Each Measurement Point' ...
              ' should be on a New Line']);
    else
        [filename,pathname] = ...
            uigetfile({'*.txt;*.dat;*.csv', ...
            'Text Files (*.txt) (*.dat) (*.csv)'} ...
            ,['Select Location List File - Each Measurement Point' ...
              ' should be on a New Line'],handles.currentDir);
    end

    if isequal(filename,0)
        disp('User selected Cancel')
        %Quit the gui
        guidata(hObject, handles);
        CloseFcn(hObject,eventdata,handles);
        return
    end

    disp(['User selected ', fullfile(pathname, filename)])

    % Open File
    FileID = fopen([pathname filename]);

    % locations is a local variable that holds location data in this
    % function
    locations = textscan(FileID,'%s','delimiter','\n');

    % append to list of reference points and convert to string array
    locations = [handles.AtlasLandmarkNames; locations{1,1}];

    % Close file
    fclose(FileID);

    % Save locations variable to be loaded next time
    disp(['Saving locations to: ', saved_location_names_loc]);
    save(saved_location_names_loc,'locations');
end

%error test the first serial port functions...
try
    %------------------------SERIAL CALLBACK SETUP---------------------
    %setup callback function to run when the polhemus system sends the
    %number of bytes assosciated with one or two sensors. NB: the stated
    %number of bytes is generally position data.
    if(handles.sensors == 1)
        handles.serial.BytesAvailableFcnCount = 48; % 48 bytes for 1 sensor
    else
        handles.serial.BytesAvailableFcnCount = 96; % 96 bytes for 2 sensors
    end
    handles.serial.BytesAvailableFcnMode = 'byte';
    handles.serial.BytesAvailableFcn = {@ReadCoordsCallback,handles};

    %--------------------------OPEN SERIAL PORT------------------------

    fopen(handles.serial);

    %-------------------STYLUS POSITION MARKING SETUP------------------
    %The following ascii string causes the stylus button to send current
    %coords. Note: 13 is the ascii code for the required newline character
    string2write = ['L1,1' 13];
    fwrite(handles.serial,string2write);
    %fwrite is used because fprintf sometimes adds extra newline characters

    %set to output in cms
    string2write = ['U1' 13];
    fwrite(handles.serial,string2write);


    %-----------------Display initial point to find on GUI-------------

    set(handles.infobox,'string',locations(1,1));

    % display locations on table in gui
    set(handles.coords_table,'Data',locations);

%catch exception if error occurs
catch serialException
    disp('COM PORT ERROR OCCURRED: Check COM Connection.')
    %run close function to close gui and delete serial port objects if
    %error occurs.
    CloseFcn(hObject,eventdata,handles);
    error(message('MATLAB:serial:fopen:opfailed', serialException.message))
end

%-------------set graph axes labels and properties---------------------

xlabel(handles.coord_plot,'X');
ylabel(handles.coord_plot,'Y');
zlabel(handles.coord_plot,'Z');


%--------------------LOAD EXPECTED COORDS---------------------------------
if ~isdeployed
    saved_expected_coords_loc = fullfile(pwd, 'saved_expected_coords.mat');
else
    saved_expected_coords_loc = fullfile(ctfroot, 'DIGIGUI', 'saved_expected_coords.mat');
end
if handles.save_expected_coords
    try
        disp(['Looking for expected_coords at: ', saved_expected_coords_loc]);
        disp(['Looking for expected_coords_tolerance at: ', saved_expected_coords_loc]);
        load(saved_expected_coords_loc);
        handles = load_expected_coords(handles, expected_coords, expected_coords_tolerance);
        handles.expected_coords = expected_coords;
        handles.expected_coords_tolerance = expected_coords_tolerance;
    catch
        uiwait(warndlg('Could not load the expected coordinates and tolerance from saved_expected_coords.mat.',...
            'Expected Coordinates Warning','modal'));
        if isfield(handles, 'expected_coords')
            handles = rmfield(handles, 'expected_coords');
        end
        if isfield(handles, 'expected_coords_tolerance')
            handles = rmfield(handles, 'expected_coords_tolerance');
        end
        handles = remove_expected_coords(handles);
    end
else
    % remove the saved file to avoid confusion
    if exist(saved_expected_coords_loc, 'file')==2
        delete(saved_expected_coords_loc);
    end
end

% Get the checkmark icon for use when checking if coordinates match
handles.checkmark_icon = imread('checkmark.tif');

% Update handles structure
guidata(hObject, handles);




% --- Import data from the given atlas directory. The directory must contain a
% file called mesh.mat which contains a struct called `mesh' within which there
% is mesh information (see the default example) and a file called
% atlas_landmarks.mat within which there is an array of landmarks called `pts'
% and a cell list of strings called `names'.
function [landmarks, landmark_names, mesh] = load_atlas_data(atlas_dir)

disp(['Loading atlas_landmarks.mat from: ', fullfile(atlas_dir, 'atlas_landmarks.mat')])
landmarks = load(fullfile(atlas_dir, 'atlas_landmarks.mat'));
landmark_names = landmarks.names;
landmarks = landmarks.pts;
disp(['Loading mesh.mat from: ', fullfile(atlas_dir, 'mesh.mat')]);
mesh = load(fullfile(atlas_dir, 'mesh.mat'));
mesh = mesh.mesh;




% --- Check atlas data matches requried specification
function [good] = good_atlas_data(landmarks, landmark_names, mesh)

good = false;
if size(landmarks, 1) < 4
    errordlg('Too few landmark points in atlas data!', 'atlas_landmarks.mat');
elseif size(landmark_names, 1) ~= size(landmarks, 1)
    errordlg('Too few/many landmark names!', 'atlas_landmarks.mat');
elseif ~all([mesh.nnode, 3] == size(mesh.node))
    errordlg('Node matrix should be nnode by 3 matrix!', 'mesh.mat');
elseif 3 ~= size(mesh.face, 2)
    errordlg('Face matrix should be an m by 3 matrix!', 'mesh.mat');
elseif any(max(mesh.face)) > mesh.nnode
    errordlg('Face matrix contains indices which exceed the dimensions of the node matrix!', 'mesh.mat');
else
    good = true;
end




% --- Executes on button press in HeadAlign.
function HeadAlign_Callback(hObject, eventdata, handles)
% hObject    handle to HeadAlign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.landmark_measurements is dynamically grown so we can always
% do an alignment when it is fully populated
if(size(handles.landmark_measurements, 1) == size(handles.AtlasLandmarks, 1))

    % extract the locations
    locations = get(handles.coords_table,'Data');
    measurements = cell2mat(locations(:, 2:4));

    %get transformation matrix to new coord system
    [TransformMatrix,TransformVector] = CoordinateTransform(handles.landmark_measurements);

    % reset list of points to just show locations to find so transformed
    % points can be plotted


    hold on

    for k = 1:size(measurements, 1)
        if isfield(handles, 'TransformMatrix')
            % undo old transform
            measurements(k,:) = measurements(k,:) * handles.TransformMatrix;
            measurements(k,:) = measurements(k,:) - handles.TransformVector;
        end

        %transform cardinal points
        measurements(k,:) = measurements(k,:) + TransformVector;
        measurements(k,:) = measurements(k,:)*TransformMatrix';

        %remove old point from graph
        delete(handles.pointhandle(k));

        %replot point
        if k <= size(handles.AtlasLandmarks, 1)
            % plot as landmark
            handles.pointhandle(k) = plot3(measurements(k,1), ...
                                           measurements(k,2), ...
                                           measurements(k,3), ...
                                           'm.', 'MarkerSize', 30, ...
                                           'Parent' , handles.coord_plot);
        else
            % plot as measurement
            handles.pointhandle(k) = plot3(measurements(k,1), ...
                                           measurements(k,2), ...
                                           measurements(k,3), ...
                                           'b.', 'MarkerSize', 30, ...
                                           'Parent' , handles.coord_plot);
        end

        %replot axes...
        axis(handles.coord_plot,'equal');

        %update newly transformed point coords (converting back
        %to a cell array first)
        locations(k,2:4) = num2cell(measurements(k,1:3));

    end

    hold off


    % Show newly transformed cardinal point coords on table
    set(handles.coords_table,'Data',locations);

    transformed_landmarks = measurements(1:size(handles.AtlasLandmarks, 1), :);

    %find matrix (A) and vector (B) needed to map head to cardinal points
    %with affine transformation
    [A,B] = affinemap(handles.AtlasLandmarks, transformed_landmarks);

    mesh_trans = handles.mesh;
    mesh_trans.node = affine_trans_RJC(handles.mesh.node,A,B);

    % Remove any existing head plot
    if isfield(handles, 'headplot')
        delete(handles.headplot)
    end

    %Then plot the transformed mesh as visual reference for further points...
    %note: this plots the head model
    hold on
    handles.headplot = trisurf(mesh_trans.face, mesh_trans.node(:,1), ...
                               mesh_trans.node(:,2),mesh_trans.node(:,3), ...
                               'FaceAlpha',0.6, ...
                               'FaceColor',[239/255 208/255 207/255], ... (skin tone rgb vals)
                               'EdgeColor','none', ...
                               'Parent',handles.coord_plot);

    % set lighting of head
    if ~isfield(handles, 'plotlight')
        handles.plotlight = light;
    end
    set(handles.headplot, 'FaceLighting', 'gouraud');

    axis equal;
    hold off;

    %save tranformation
    handles.TransformMatrix = TransformMatrix;
    handles.TransformVector = TransformVector;

    % Update handles structure
    guidata(hObject, handles);

end

function save_settings(handles)

if ~isdeployed
    settings_loc = fullfile(pwd, 'settings.mat');
else
    settings_loc = fullfile(ctfroot, 'DIGIGUI', 'settings.mat');
end

% save settings
if isfield(handles, 'atlas_dir')
    disp(['Saving atlas_dir to ', settings_loc]);
    atlas_dir = handles.atlas_dir;
    save(settings_loc, 'atlas_dir');
end
if isfield(handles, 'error_distance')
    disp(['Saving error_distance to ', settings_loc]);
    error_distance = handles.error_distance;
    save(settings_loc, 'error_distance', '-append');
end
if isfield(handles, 'double_tap_error_enabled')
    disp(['Saving double_tap_error_enabled to ', settings_loc]);
    double_tap_error_enabled = handles.double_tap_error_enabled;
    save(settings_loc, 'double_tap_error_enabled', '-append');
end
if isfield(handles, 'save_expected_coords')
    disp(['Saving save_expected_coords to ', settings_loc]);
    save_expected_coords = handles.save_expected_coords;
    save(settings_loc, 'save_expected_coords', '-append');
end
if isfield(handles, 'only_measure_with_expected_coords')
    disp(['Saving only_measure_with_expected_coords to ', settings_loc]);
    only_measure_with_expected_coords = handles.only_measure_with_expected_coords;
    save(settings_loc, 'only_measure_with_expected_coords', '-append');
end
disp(['Saving insert_by_location_name to ', settings_loc]);
insert_by_location_name = handles.menu_options_insert_by_location_name.Checked
save(settings_loc, 'insert_by_location_name', '-append');
if isfield(handles, 'previousStartNumber')
    disp(['Saving previousStartNumber to ', settings_loc]);
    previousStartNumber = handles.previousStartNumber;
    save(settings_loc, 'previousStartNumber', '-append');
end
if isfield(handles, 'previousEndNumber')
    disp(['Saving previousEndNumber to ', settings_loc]);
    previousEndNumber = handles.previousEndNumber;
    save(settings_loc, 'previousEndNumber', '-append');
end
if isfield(handles, 'previousStartLetter')
    disp(['Saving previousStartLetter to ', settings_loc]);
    previousStartLetter = handles.previousStartLetter;
    save(settings_loc, 'previousStartLetter', '-append');
end
if isfield(handles, 'previousEndLetter')
    disp(['Saving previousEndLetter to ', settings_loc]);
    previousEndLetter = handles.previousEndLetter;
    save(settings_loc, 'previousEndLetter', '-append');
end


function save_locations(handles)

choice = questdlg('Do you want to save the locations list for loading next time or reset the list to what was first imported upon starting this session?', ...
    'Locations List', ...
    'Save Locations List','Reset Locations List','Reset Locations List');
if strcmp(choice, 'Save Locations List')
    % Save locations variable to be loaded next time
    if ~isdeployed
        saved_location_names_loc = fullfile(pwd,'savedLocationNames.mat');
    else
        saved_location_names_loc = fullfile(ctfroot, 'DIGIGUI', 'savedLocationNames.mat');
    end
    locations = handles.coords_table.Data(:, 1);
    disp(['Saving locations to: ', saved_location_names_loc]);
    save(saved_location_names_loc,'locations');
end


function handles = close_serial_port(handles)

%close port only if not closed
if(isfield(handles,'COMport'))
    if(handles.COMport ~= 0)
        if( ~ strcmp(handles.serial.status, 'closed') )
            fclose(handles.serial);
        end
        delete(handles.serial);
    end
end

return


function CloseFcn(source,event,handles)
%my user-defined close request function

try
    handles = guidata(handles.figure1);
catch
    handles = struct([]);
end

save_settings(handles);
save_locations(handles);
close_serial_port(handles);

delete(gcf);


function ReadCoordsCallback(s,BytesAvailable,handles)

% Update handles structure to most current version
handles = guidata(handles.figure1);

%read the data on the serial port that triggered the callback
data_str=fgetl(s);

%read a second line if there are two sensors
if(handles.sensors == 2)
    data_str(2,:) = fgetl(s);
end

if isfield(handles, 'disable_measurements') && handles.disable_measurements
    return
end
% Don't measure if we have a double tap error open
if isfield(handles, 'doubleTapErrorFigure') && isvalid(handles.doubleTapErrorFigure)
    return
end
if (...
        isfield(handles, 'only_measure_with_expected_coords') && ...
        handles.only_measure_with_expected_coords && ...
        ~isfield(handles, 'expected_coords') ...
)
    msg = 'Measurements cannot be made until a list of expected coordinates has been imported! This error can be disabled in the settings.';
    errordlg(msg, 'Missing Expected Coordinates', 'modal');
    return
end


data_num=str2num(data_str);

% Format of data obtained for the current settings
% 1 2 3 4 5 6 7
% 1 Detector Number (should be 1 for 1 stylus)
% 2 X position in cms
% 3 Y position in cms
% 4 Z position in cms
% 5 Azimuth of stylus in degrees
% 6 Elevation of stylus in degrees
% 7 Roll of stylus degrees

% extract coords
Coords = data_num(:,2:4);

% if there are 2 sensors do vector subtraction to get position of
% stylus sensor relative to second sensor
if(handles.sensors == 2)
    Coords = Coords(1,:) - Coords(2,:);
end

% Extract previous data from table
data = get(handles.coords_table,'Data');

%increment the point count to the next unmeasured point before measurement
previous_point_count = handles.point_count;
at_unmeasured_point = false;
while ~at_unmeasured_point
    try
        handles.point_count = handles.point_count + 1;
        next_x_coord = data(handles.point_count, 2);
        if isempty(next_x_coord{1})
            at_unmeasured_point = true;
        end
    catch
        % will get 'Index exceeds matrix dimensions' error for very first
        % measurement since we have no cells in column 2 yet
        at_unmeasured_point = true;
    end
end
point_count_increment = handles.point_count - previous_point_count;

% Save original coordinates of landmark positions if measuring those
if(handles.point_count <= size(handles.AtlasLandmarks, 1))
    handles.landmark_measurements(handles.point_count, :) = Coords;
end
% disable head alignment butten until we have all landmark positions.
% Note that landmark_measurements is a dynamically resized array so
% the below check works!
if ~isfield(handles, 'landmark_measurements') || size(handles.landmark_measurements, 1) < size(handles.AtlasLandmarks, 1)
    set(handles.HeadAlign,'Enable','off');
else
    set(handles.HeadAlign,'Enable','on');
end
% Do coord transform on points measured after alignment - i.e. if
% there is a head plotted.
if isfield(handles, 'headplot')
    Coords = Coords + handles.TransformVector;
    Coords = Coords*handles.TransformMatrix';
end

location_list_expanded = false;
% Check if table is currently full - if it is then adding a new point
% will expand the table...
if(handles.point_count > size(data,1))
    % ... so update the bool that tracks if location names have been
    % edited. When the user saves their data they will therefore be
    % prompted to save the locations list too.
    handles.editedLocationsList = true;
    location_list_expanded = true;
end

% Update table with newly measured x y and z values
data(handles.point_count,2:4) = num2cell(Coords);
set(handles.coords_table,'Data',data);

% Double tap warning...
if handles.double_tap_error_enabled && handles.point_count > 1
    % extract the previous measured point for comparison
    last_point = cell2mat(data(handles.point_count-point_count_increment, 2:4));
    % If our last point is empty then we can't do a comparison. This
    % will happen when you select a different row to measure, in which
    % case a double tap isn't going to happen
    if length(last_point) == 3
        distance = norm(Coords - last_point);
        if distance < handles.error_distance
            last_point = data(handles.point_count-point_count_increment, 1);
            this_point = data(handles.point_count, 1);
            msg = sprintf('%s measurement was only %0.2g cm from %s measurement!\nCurrent error distance is %0.2g cm. This can be changed or disabled in the options.', this_point{1}, distance, last_point{1}, handles.error_distance);
            handles.doubleTapErrorFigure = errordlg(msg, 'Double tap error', 'modal');
            if location_list_expanded
                % delete row
                data(handles.point_count, :) = [];
            else
                % remove data, restore measurement status and distance
                data(handles.point_count,2:4) = {[], [], []};
                if size(data, 2) > 4 && isfield(handles, 'expected_coords')
                    measurement_status = data(handles.point_count, 8);
                    measurement_status = measurement_status{1};
                    if ~isempty(measurement_status)
                        data(handles.point_count, 8) = {'Unmeasured'};
                    end
                    data(handles.point_count, 10) = {[]};
                end
            end
            set(handles.coords_table,'Data',data);
            % reset point count, update guidata and exit
            handles.point_count = previous_point_count;
            guidata(handles.figure1,handles);
            return
        end
    end
end

% Remove any open expected or unexpected measurement figures
if isfield(handles, 'expectedMeasurementFigure') && isvalid(handles.expectedMeasurementFigure)
    close(handles.expectedMeasurementFigure)
end
if isfield(handles, 'unexpectedMeasurementWarnFigure') && isvalid(handles.unexpectedMeasurementWarnFigure)
    close(handles.unexpectedMeasurementWarnFigure)
end
% Check against expected coordinates if any are specified
if isfield(handles, 'expected_coords')
    rowmatch = find(strcmp(data{handles.point_count, 1}, handles.expected_coords.Properties.RowNames));
    if ~isempty(rowmatch)
        distance = norm(Coords - handles.expected_coords(rowmatch, :).Variables);
        this_point = data(handles.point_count, 1);
        if distance > handles.expected_coords_tolerance
            msg = sprintf('%s measurement is %0.2g cm from the expected location of (%0.2g, %0.2g, %0.2g) cm!\nCurrent tolerance is %0.2g cm. Tolerance is set in the expected coordinates file.', this_point{1}, distance, handles.expected_coords(rowmatch, 1).Variables, handles.expected_coords(rowmatch, 2).Variables, handles.expected_coords(rowmatch, 3).Variables, handles.expected_coords_tolerance);
            handles.unexpectedMeasurementWarnFigure = warndlg(msg, 'Unexpected Measurement!');
            % Set measurement status...
            data(handles.point_count, 8) = {'No'};
        else
            msg = sprintf('%s measurement is within tolerance of (%0.2g, %0.2g, %0.2g) cm.\nCurrent tolerance is %0.2g cm. Tolerance is set in the expected coordinates file.', this_point{1}, handles.expected_coords(rowmatch, 1).Variables, handles.expected_coords(rowmatch, 2).Variables, handles.expected_coords(rowmatch, 3).Variables, handles.expected_coords_tolerance);
            handles.expectedMeasurementFigure = msgbox(msg, 'Expected Measurement Success!', 'custom', handles.checkmark_icon);
            % Set measurement status...
            data(handles.point_count, 8) = {'Yes'};
        end
        % Display distance
        data(handles.point_count, 10) = {distance};
        set(handles.coords_table,'Data',data);
    end
end

% update point to look for (unless at end of list as given by the
% length of data - ie the number of headpoints)
if( handles.point_count < size(data,1) )
    % search for next unmeasured point
    point_count_increment = 0;
    at_unmeasured_point = false;
    while ~at_unmeasured_point
        try
            point_count_increment = point_count_increment + 1;
            next_x_coord = data(handles.point_count+point_count_increment, 2);
            if isempty(next_x_coord{1})
                at_unmeasured_point = true;
            end
        catch
            % will get 'Index exceeds matrix dimensions' error for very first
            % measurement since we have no cells in column 2 yet
            at_unmeasured_point = true;
        end
    end
    if handles.point_count+point_count_increment <= size(data,1)
        set(handles.infobox,'string', data(handles.point_count+point_count_increment,1));
    else
        set(handles.infobox,'string','End of locations list reached');
    end
else
    set(handles.infobox,'string','End of locations list reached');
end


% Remove point from graph if present...
if isfield(handles, 'pointhandle')
    if handles.point_count <= size(handles.pointhandle, 2)
        delete(handles.pointhandle(handles.point_count));
        % and replot axes.
        axis(handles.coord_plot,'equal');
    end
end

%add the measured point to the 3d graph
hold(handles.coord_plot,'on');
%save the handle of the point so it can be removed later...
if(handles.point_count <= size(handles.AtlasLandmarks, 1))
    handles.pointhandle(handles.point_count) = plot3(Coords(1), ...
                                        Coords(2),Coords(3), ...
                                        'm.', 'MarkerSize', 30, ...
                                        'Parent' , handles.coord_plot);
else %Note: above marker points are plotted differently
    handles.pointhandle(handles.point_count) = plot3(Coords(1), ...
                                        Coords(2),Coords(3), ...
                                       'b.', 'MarkerSize', 30, ...
                                       'Parent',handles.coord_plot);
end
hold(handles.coord_plot,'off');
%replot axes...
axis(handles.coord_plot,'equal');

% save the point_count so it can be undone
handles.point_count_history(end+1) = handles.point_count;

% Update handles structure
guidata(handles.figure1,handles);


% --- Executes on button press in remove_last_pt.
function remove_last_pt_Callback(hObject, eventdata, handles)
% hObject    handle to remove_last_pt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%don't delete points if alignment already done or at first point
if ~isempty(handles.point_count_history)

    data = get(handles.coords_table,'Data');

    last_point_count = handles.point_count_history(end);
    handles.point_count_history(end) = [];

    % remove data, restore measurement status and distance
    data(last_point_count,2:4) = {[], [], []};
    if size(data, 2) > 4 && isfield(handles, 'expected_coords')
        measurement_status = data(last_point_count, 8);
        measurement_status = measurement_status{1};
        if ~isempty(measurement_status)
            data(last_point_count, 8) = {'Unmeasured'};
        end
        data(last_point_count, 10) = {[]};
    end
    set(handles.coords_table,'Data',data);

    % Remove point from graph...
    delete(handles.pointhandle(last_point_count));
    % and replot axes.
    axis(handles.coord_plot,'equal');

    % Reset point_count so next measurement is of the point which
    % has just been deleted
    handles.point_count = last_point_count-1;

    data = get(handles.coords_table,'Data');

    % update next point to look for string
    set(handles.infobox,'string', data(handles.point_count+1,1));

    % Disable align if now not enough points
    if ~all(ismember(1:size(handles.AtlasLandmarks, 1), handles.point_count_history))
        set(handles.HeadAlign,'Enable','off');
    end

    % Update handles structure
    guidata(handles.figure1,handles);

end



% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% disable measurements
handles.disable_measurements = true;
guidata(hObject,handles);

% Find interface objects that are set to 'on' i.e. enabled...
InterfaceObj=findobj(handles.figure1,'Enable','on');
% ... and turn them off.
set(InterfaceObj,'Enable','off');

% Open a "Save As..." Dialogue with different saving options as shown.
% The filterIndex gives the index (1, 2 or 3) of the chosen save type.
if ~isdeployed
    [fileName,pathName,filterIndex] = ...
        uiputfile({'*.csv;*.dat;*.txt', ...
        'Comma-delimited text files (*.csv) (*.dat) (*.txt)'; ...
        ...
        '*.mat', ...
        'MAT-file (*.mat)'; ...
        ...
        '*.xls;*.xlsb;*.xlsm;*.xlsx', ...
        'Excel� spreadsheet files (*.xls) (*.xlsb) (*.xlsm) (*.xlsx)'; ...
        },'Save As...');
else
    [fileName,pathName,filterIndex] = ...
        uiputfile({'*.csv;*.dat;*.txt', ...
        'Comma-delimited text files (*.csv) (*.dat) (*.txt)'; ...
        ...
        '*.mat', ...
        'MAT-file (*.mat)'; ...
        ...
        '*.xls;*.xlsb;*.xlsm;*.xlsx', ...
        'Excel� spreadsheet files (*.xls) (*.xlsb) (*.xlsm) (*.xlsx)'; ...
        },'Save As...',handles.userDir);
end

% Re-enable the interface objects.
set(InterfaceObj,'Enable','on');

% re-enable measurements
handles.disable_measurements = false;
guidata(hObject,handles);

if(filterIndex ~= 0) % if == 0 then user selected "cancel" in "Save As"

    data = get(handles.coords_table,'Data');

    % check data cell array has same number of columns as there are column
    % names.
    if(size(data,2) < length(get(handles.coords_table,'ColumnName')))

        % dont save table if not enough data is available
        errordlg('Cannot save without recorded location data.', ...
            'Save Error','modal');

        % exit function here
        guidata(hObject,handles);
        return;
    end

    % If the chosen save type is .mat then use a standard matlab save command
    if(filterIndex == 2)
        disp(['Data saving to ' pathName fileName]);
        dataOutput = get(handles.coords_table,'Data');
        dataOutputHeadings = handles.coords_table.ColumnName;
        save([pathName fileName],'dataOutput', 'dataOutputHeadings');
        disp('Data is stored in cell array "dataOutput" with headings in "dataOutputHeadings"');

    % Otherwise create a table from the cell array and output that to file.
    else
        % find any empty cells in Locations data
        emptyLocationNames = cellfun('isempty',data(:,1));
        buttonPressed = 'Yes';

        if(any(emptyLocationNames))
            % Warn the user if there are any location names missing...
            buttonPressed = questdlg({'Some location names are unspecified.';
                                      'Missing location names will be replaced by the symbol "-".';...
                                      'Would you like to continue?'},...
                                      'Warning','Yes','No','Yes');
        end

        %Only save data if user presses Yes or Yes has been set previously.
        if(strcmp(buttonPressed,'Yes'))

            disp(['Data saving to ' pathName fileName]);

            %Mark empty location names as '-'
            data(emptyLocationNames,1) = {'-'};

            % replace spaces with underscores in column names
            columnnames = regexprep(handles.coords_table.ColumnName, ' +', '_');

            tableToOutput = cell2table(data,'VariableNames', columnnames);
            % Note that writetable changes its output depending on the fileName
            % type.
            writetable(tableToOutput,[pathName fileName]);
        end
    end

    if(handles.editedLocationsList) % true if edited
        % check if user wants to save locations list too
        button = 'No';
        button = questdlg({['The locations list appears to have been'...
            ' edited or added to since it was last imported.'];...
            ['Would you also like to export your current '...
            'locations list?']},...
            'Export?','Yes','No','Yes');
        if(strcmp(button,'Yes'))
            % Call the export headpoints callback
            handles = guidata(handles.figure1);
            ExportHeadpoints_Callback(hObject, eventdata, handles);
        end
    end

end
guidata(hObject,handles);


% --- Executes when selected cell(s) is changed in coords_table.
function coords_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to coords_table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

if(isempty(eventdata.Indices))
    % delete 'selectedRow' field of 'handles' if the callback is triggered
    % by deselection (eg by the removal or addition of a set of coordinates)
    handles = rmfield(handles,'selectedRow');
else
    % extract the row(s) from where the user clicked on the table.
    handles.selectedRow = eventdata.Indices(:,1);
end
guidata(hObject,handles);



% --- Executes on button press in InsertRowPushbutton.
function InsertRowPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to InsertRowPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
insert_rows_below(hObject, eventdata, handles, 1)



% --- Executes on button press in DeleteRowPushbutton.
function DeleteRowPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteRowPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = get(handles.coords_table,'Data');

% See if the selectedRow variable exists within the handles struct
% (doesn't if no selection performed before clicking or cell has been
% deselected)
if(isfield(handles,'selectedRow'))
    if(handles.selectedRow(1) <= size(handles.AtlasLandmarks, 1))
        errordlg('Cannot insert or delete Atlas Points','Edit Error','modal');
    else
        % delete selected rows...
        data(handles.selectedRow,:) = [];

        % Locations list has now been edited so change bool.
        handles.editedLocationsList = true;

        % make sure selected rows are in descending order so we can remove
        % them from the point count history (we decrement as necessary)
        descendingSelectedRows = sort(handles.selectedRow, 'descend');

        numToDecrement = 0;
        for i=1:length(descendingSelectedRows)
            selectedRow = descendingSelectedRows(i);
            [isin, idx] = ismember(selectedRow, handles.point_count_history);
            if any(isin)
                % remove point and decrement point numbers above
                handles.point_count_history(idx) = [];
                to_decrement = handles.point_count_history > selectedRow;
                handles.point_count_history(to_decrement) = handles.point_count_history(to_decrement) - 1;
                % Remove point from graph
                delete(handles.pointhandle(selectedRow));
                % will need to change point count...
                numToDecrement = numToDecrement + 1;
            end
        end
        if numToDecrement > 0
            % replot axes
            axis(handles.coord_plot,'equal');
        end
        handles.point_count = handles.point_count - numToDecrement;
    end
else
    % Tell user to select a row before inserting
    errordlg('Please select a row to delete.','Delete Error','modal');
end
% save the newly changed data to the table on the gui
set(handles.coords_table,'Data',data);

guidata(hObject,handles);


% --- Executes on button press in ImportHeadpoints.
function ImportHeadpoints_Callback(hObject, eventdata, handles)
% hObject    handle to ImportHeadpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% disable measurements
handles.disable_measurements = true;
guidata(hObject,handles);

% Find interface objects that are set to 'on' i.e. enabled...
InterfaceObj=findobj(handles.figure1,'Enable','on');
% ... and turn them off.
set(InterfaceObj,'Enable','off');

%--------------------HEADPOINTS TO DIGITISE INPUT-----------------------
if ~isdeployed
    [filename,pathname] = ...
        uigetfile({'*.txt;*.dat;*.csv', ...
        'Text Files (*.txt) (*.dat) (*.csv)'} ...
        ,['Select Location List File - Each Measurement Point Should be'...
        ' on a New Line']);
else
    [filename,pathname] = ...
        uigetfile({'*.txt;*.dat;*.csv', ...
        'Text Files (*.txt) (*.dat) (*.csv)'} ...
        ,['Select Location List File - Each Measurement Point Should be'...
        ' on a New Line'],handles.userDir);
end
% Re-enable the interface objects.
set(InterfaceObj,'Enable','on');

% re-enable measurements
handles.disable_measurements = false;
guidata(hObject,handles);

% user selected cancel...
if isequal(filename,0)
    return
end

% Warn user that this will reset all currently gathered data if any has
% been collected.
if(handles.point_count > 0)

    button = 'No';

    button = questdlg({'Warning! Any existing data will be lost.';...
        'Do you wish to continue?'},'Data Warning','Yes','No','No');

    % user selected cancel...
    if strcmp(button,'No')
        return
    end

end
% Warn user that this will reset any loaded expected coordinates
if isfield(handles, 'expected_coords') || isfield(handles, 'expected_coords_tolerance')

    button = questdlg({'Warning! Any loaded expected coordinates will be lost.';...
        'Do you wish to continue?'},'Expected Coordinates Warning','Yes','No','No');

    % user selected cancel...
    if strcmp(button,'No')
        return
    end
end

if isfield(handles, 'expected_coords')
    handles = rmfield(handles, 'expected_coords');
end
if isfield(handles, 'expected_coords_tolerance')
    handles = rmfield(handles, 'expected_coords_tolerance');
end
handles = remove_expected_coords(handles);

disp(['User selected ', fullfile(pathname, filename)])

% Open File
FileID = fopen([pathname filename]);

% locations is a local variable that holds location data in this
% function
locations = textscan(FileID,'%s','delimiter','\n');

% append to list of reference points and convert to string array
locations = [handles.AtlasLandmarkNames; locations{1,1}];

% Close file
fclose(FileID);

% Reset points counter
handles.point_count = 0;

% Display initial point to find on GUI
set(handles.infobox,'string',locations(1,1));

% display locations on table in gui
set(handles.coords_table,'Data',locations);

% if head align button has been enabled set to disabled.
if(strcmp(get(handles.HeadAlign,'Enable'),'on'))
    set(handles.HeadAlign,'Enable','off');
end

% clear previous measurements and headmap from plot...
cla(handles.coord_plot);
% and replot axes.
axis(handles.coord_plot,'equal');

% Locations list is now unedited (since it's just been imported) so reset
% both bools that deal with whether bits of the location list are edited.
handles.editedLocationsList = false;
handles.editedAtlasPoints = false;

guidata(hObject,handles);


% --- Executes on button press in ExportHeadpoints.
function ExportHeadpoints_Callback(hObject, eventdata, handles)
% hObject    handle to ExportHeadpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Find interface objects that are set to 'on' i.e. enabled...

% disable measurements
handles.disable_measurements = true;
guidata(hObject,handles);

% Find interface objects that are set to 'on' i.e. enabled...
InterfaceObj=findobj(handles.figure1,'Enable','on');
% ... and turn them off.
set(InterfaceObj,'Enable','off');

% Check that user is happy to continue
button = questdlg('Warning! This only exports location names. To save location coordinates use "Save As...". Do you wish to continue?' ...
    ,'Export Warning','Yes','No','No');
if strcmp(button, 'No')
    % Re-enable the interface objects.
    set(InterfaceObj,'Enable','on');
    % re-enable measurements
    handles.disable_measurements = false;
    guidata(hObject,handles);
    return
end

% Set name to previous name prior to editing if user selects "no"
if(strcmp(button,'No'))
    data = get(handles.coords_table,'Data');
    data{selectedRow,1} = PreviousData;
    set(handles.coords_table,'Data',data);
    NewData = PreviousData;
else
    % Atlas points have been edited so update bool.
    handles.editedAtlasPoints = true;
end

% Display warning dialogue before uiputfile if the atlas points have been
% editied...
if(handles.editedAtlasPoints)
    uiwait(warndlg({['Note: Atlas points are NOT included in location'...
        ' list files.'];...
        ['Any atlas point renaming will NOT be reflected the'...
        ' exported file.']}...
    ,'Export Warning','modal'));
end

% Open an "Export" Dialogue
if ~isdeployed
    [fileName,pathName,filterIndex] = ...
        uiputfile({'*.txt;*.dat;*.csv', ...
        'Text Files (*.txt) (*.dat) (*.csv)'} ...
        ,'Export Location List File ...');
else
    [fileName,pathName,filterIndex] = ...
        uiputfile({'*.txt;*.dat;*.csv', ...
        'Text Files (*.txt) (*.dat) (*.csv)'} ...
        ,'Export Location List File ...',handles.userDir);
end

% Re-enable the interface objects.
set(InterfaceObj,'Enable','on');

% re-enable measurements
handles.disable_measurements = false;
guidata(hObject,handles);

% Otherwise create a table from the cell array and output that to file.
if(filterIndex ~= 0) % if == 0 then user selected "cancel" in save dialogue

    data = get(handles.coords_table,'Data');

    % error if outputting only atlas points
    if(size(data,1) <= size(handles.AtlasLandmarks, 1))
        errordlg({'Cannot export locations:';...
            'Only atlas point locations have been found.';...
            'Atlas points alone cannot be exported.'},...
            'Export Error','modal');
    else

        disp(['Locations saving to ' pathName fileName]);

        fileID = fopen([pathName fileName],'wt');

        %write from after the atlas points to the last data point
        for i = size(handles.AtlasLandmarks, 1)+1:size(data,1)
            fprintf(fileID,'%s\n',data{i,1});
        end

        fclose(fileID);
        clear fileID;

    end

end

guidata(hObject,handles);


% --- Executes on button press in measureThisRowButton.
function measureThisRowButton_Callback(hObject, eventdata, handles)
% hObject    handle to measureThisRowButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% See if the selectedRow variable exists within the handles struct
% (doesn't if no selection performed before clicking or cell has been
% deselected)
if(isfield(handles,'selectedRow'))

    % Set point_count such that the first of the selected rows will be
    % measured
    handles.point_count = handles.selectedRow(1)-1;

    % Update the "Point to Get" string
    data = get(handles.coords_table,'Data');
    set(handles.infobox,'string',...
            data(handles.selectedRow(1),1))

    % delete the point data if there's data to delete
    if size(data, 2) > 1
        % Remove plotted points from graph and replot axes
        for i = 1:length(handles.selectedRow)
            row_x_coord = data(handles.selectedRow(i), 2);
            if ~isempty(row_x_coord{1})
                delete(handles.pointhandle(handles.selectedRow(i)));
                axis(handles.coord_plot,'equal');
                [isin, idx] = ismember(handles.selectedRow(i), handles.point_count_history);
                assert(all(isin));
                handles.point_count_history(idx) = [];
            end
        end
        % delete data on table
        data(handles.selectedRow,2:4) = cell(length(handles.selectedRow), 3);
        if size(data, 2) > 4
            % Reset measurement status and distance if we have expected measurements
            for i = 1:length(handles.selectedRow)
                selectedRow = handles.selectedRow(i);
                measurement_status = data(selectedRow, 8);
                measurement_status = measurement_status{1};
                if ~isempty(measurement_status)
                    data(selectedRow, 8) = {'Unmeasured'};
                end
                data(selectedRow, 10) = {[]};
            end
        end
        set(handles.coords_table,'Data',data);
    end

else
    % No point selected
    errordlg('Please select rows to delete and gather location data.',...
        'Selection Error','modal');
end
guidata(hObject,handles);


% --- Executes when entered data in editable cell(s) in coords_table.
function coords_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to coords_table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

% Extract selected row
selectedRow = eventdata.Indices(:,1);

% Extract previous and new data (the eventdata fields are read only hence
% writing them to a new variable).
NewData = eventdata.NewData;
PreviousData = eventdata.PreviousData;

% warn when editing of rows less than number of atlas points
if(selectedRow <= size(handles.AtlasLandmarks, 1))

    % Check that user is happy to continue
    button = 'Yes';
    button = questdlg({['Warning! About to rename Atlas Point "' ...
        PreviousData, '" to "', NewData, '".']; ...
        '';...
        'Any changes will NOT be reflected in Exported Locations files';...
        '(exported using the "Export Locations..." button) but WILL be';...
        ['reflected in saved data files (saved using the "Save Data'...
        ' As..." button.)'];...
        '';...
        'Do you wish to continue?'} ...
        ,'Rename Warning','Yes','No','No');

    % Set name to previous name prior to editing if user selects "no"
    if(strcmp(button,'No'))
        data = get(handles.coords_table,'Data');
        data{selectedRow,1} = PreviousData;
        set(handles.coords_table,'Data',data);
        NewData = PreviousData;
    else
        % Atlas points have been edited so update bool.
        handles.editedAtlasPoints = true;
    end

end
% if the previous and the new data are different set editedLocationsList to
% true.
if(~strcmp(PreviousData,NewData))
    handles.editedLocationsList = true;
end

guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_file_import_atlas_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_import_atlas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% disable measurements
handles.disable_measurements = true;
guidata(hObject,handles);

% Find interface objects that are set to 'on' i.e. enabled...
InterfaceObj=findobj(handles.figure1,'Enable','on');
% ... and turn them off.
set(InterfaceObj,'Enable','off');

%--------------------DIRECTORY INPUT-----------------------
atlas_dir = uigetdir(handles.userDir, ['Select atlas directory']);

% Re-enable the interface objects.
set(InterfaceObj,'Enable','on');

% re-enable measurements
handles.disable_measurements = false;
guidata(hObject,handles);

% user selected cancel...
if isequal(atlas_dir,0)
    return
end

% Warn user that this will reset all currently gathered data if any has
% been collected.
if(handles.point_count > 0)

    button = 'No';

    button = questdlg({'Warning! Any existing data will be lost.';...
        'Do you wish to continue?'},'Data Warning','Yes','No','No');

    % user selected cancel...
    if strcmp(button,'No')
        return
    end

end

% Check atlas data is good
noatlas = true;
while noatlas
    try
        %load other data needed for headpoint plotting - these are the required
        %names after getting atlas_dir from settings.mat
        [landmarks, landmark_names, mesh] = load_atlas_data(atlas_dir);
        if good_atlas_data(landmarks, landmark_names, mesh)
            noatlas = false;
        else
            return
        end
    catch
        errordlg('Bad atlas directory!', atlas_dir);
        return
    end
end

% save our new mesh and landmarks
handles.AtlasLandmarks = landmarks;
handles.AtlasLandmarksNames = landmark_names;
handles.mesh = mesh;

% Remove the old atlas directory from the path and add the new one to get
% necessary functions
rmpath(handles.atlas_dir)
handles.atlas_dir = atlas_dir;
addpath(handles.atlas_dir);

% Reset points counter
handles.point_count = 0;

% The landmark names are now the nes locations which we save
locations = landmark_names;

% Display initial point to find on GUI
set(handles.infobox,'string',locations(1,1));

% display locations on table in gui
set(handles.coords_table,'Data',locations);

% if head align button has been enabled set to disabled.
if(strcmp(get(handles.HeadAlign,'Enable'),'on'))
    set(handles.HeadAlign,'Enable','off');
end

% clear previous measurements and headmap from plot...
cla(handles.coord_plot);
% and replot axes.
axis(handles.coord_plot,'equal');

guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_import_locations_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_import_locations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ImportHeadpoints_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function menu_file_export_locations_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_export_locations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ExportHeadpoints_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_options_Callback(hObject, eventdata, handles)
% hObject    handle to menu_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_options_double_tap_Callback(hObject, eventdata, handles)
% hObject    handle to menu_options_double_tap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% disable measurements
handles.disable_measurements = true;
guidata(hObject,handles);

% Find interface objects that are set to 'on' i.e. enabled...
InterfaceObj=findobj(handles.figure1,'Enable','on');
% ... and turn them off.
set(InterfaceObj,'Enable','off');

%--------------------DIALOGUE BOX-----------------------
[handles.error_distance, handles.double_tap_error_enabled] = ...
    doubletapdialog(handles.error_distance, handles.double_tap_error_enabled);

% Re-enable the interface objects.
set(InterfaceObj,'Enable','on');

% re-enable measurements
handles.disable_measurements = false;
guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_edit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_edit_undo_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
remove_last_pt_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_file_saveas_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_saveas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
save_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_file_import_expected_coordinates_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_import_expected_coordinates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% disable measurements
handles.disable_measurements = true;
guidata(hObject,handles);

% Find interface objects that are set to 'on' i.e. enabled...
InterfaceObj=findobj(handles.figure1,'Enable','on');
% ... and turn them off.
set(InterfaceObj,'Enable','off');

%--------------------HEADPOINTS TO DIGITISE INPUT-----------------------
if ~isdeployed
    [filename,pathname] = ...
        uigetfile({'*.txt;*.dat;*.csv', ...
        'Text Files (*.txt) (*.dat) (*.csv)'} ...
        ,['Select expected coordinates list file - each location and coordinates should be'...
        ' on a new line. Tolerance is taken as the X coordinate of the "Tolerance" location.']);
else
    [filename,pathname] = ...
        uigetfile({'*.txt;*.dat;*.csv', ...
        'Text Files (*.txt) (*.dat) (*.csv)'} ...
        ,['Select expected coordinates list file - each location and coordinates should be'...
        ' on a new line. Tolerance is taken as the X coordinate of the "Tolerance" location.'],handles.userDir);
end
% Re-enable the interface objects.
set(InterfaceObj,'Enable','on');

% re-enable measurements
handles.disable_measurements = false;
guidata(hObject,handles);

% user selected cancel...
if isequal(filename,0)
    return
end

disp(['User selected ', fullfile(pathname, filename)])

expected_coords = readtable([pathname filename]);

% Check file specification
if ~strcmp(expected_coords.Properties.VariableNames{1}, 'Location')
    errordlg('Import Error: First entry in top line of coordinates list file should be "Location" (note capitalisation).', 'Import Error');
    return
end
if ~strcmp(expected_coords.Properties.VariableNames{2}, 'X')
    errordlg('Import Error: Second entry in top line of coordinates list file should be "X" (note capitalisation).', 'Import Error');
    return
end
if ~strcmp(expected_coords.Properties.VariableNames{3}, 'Y')
    errordlg('Import Error: Third entry in top line of coordinates list file should be "Y" (note capitalisation).', 'Import Error');
    return
end
if ~strcmp(expected_coords.Properties.VariableNames{4}, 'Z')
    errordlg('Import Error: Fourth entry in top line of coordinates list file should be "Z" (note capitalisation).', 'Import Error');
    return
end

% Find tolerance
tolrow = find(strcmp(expected_coords.Location, 'Tolerance'));
if isempty(tolrow)
    errordlg('Tolerance not found in expected coordinates file! It should be an entry in the "Location" column.', 'Import Error');
    return
end
if length(tolrow) > 1
    errordlg('Multiple tolerance entries found in expected coordinates file! There should only be one.', 'Import Error');
    return
end
try
    expected_coords_tolerance = expected_coords.X(tolrow);
    assert(isnumeric(expected_coords_tolerance));
    assert(expected_coords_tolerance > 0);
catch
    errordlg('Could not get tolerance from expected coordinates file! It should be a nonnegative number in the "X" column next to the "Tolerance" location.', 'Import Error');
    return
end

% Delete tolerance from location names
expected_coords.Properties.RowNames = expected_coords.Location;
expected_coords.Location = [];
expected_coords('Tolerance', :) = [];

handles = load_expected_coords(handles, expected_coords, expected_coords_tolerance);

% save
handles.expected_coords = expected_coords;
handles.expected_coords_tolerance = expected_coords_tolerance;

% Save locations variable to be loaded next time
if ~isdeployed
    saved_expected_coords_loc = fullfile(pwd, 'saved_expected_coords.mat');
else
    saved_expected_coords_loc = fullfile(ctfroot, 'DIGIGUI', 'saved_expected_coords.mat');
end
disp(['Saving expected_coords to: ', saved_expected_coords_loc]);
disp(['Saving expected_coords_tolerance to: ', saved_expected_coords_loc]);
save(saved_expected_coords_loc,'expected_coords','expected_coords_tolerance');

guidata(hObject,handles);


function handles = load_expected_coords(handles, expected_coords, expected_coords_tolerance)
% handles    structure with handles and user data (see GUIDATA)
% expected_coords loaded expected coordinates table with tolerance removed

% Search for expected coordinates names in location names
data = get(handles.coords_table,'Data');
locations = data(:, 1);
locations_with_names = locations(~cellfun('isempty',locations));
[is_match, match_rows] = ismember(expected_coords.Properties.RowNames, locations_with_names);
if ~all(is_match)
    errordlg('Not all expected coordinate locations found in currently loaded coordinate locations! Check the expected coordinates file and the list of coordinates.', 'Import Error');
    return
end

% Add as new rows to data in new columns
if size(data, 2) == 1
    % no data recorded yet
    data(:, 2:4) = cell(size(data,1), 3);
end
data(:, 5:7) = cell(size(data,1), 3);
data(match_rows, 5:7) = table2cell(expected_coords);

% Add measurement status too
data(:, 8) = cell(size(data,1), 1);
data(match_rows, 8) = {'Unmeasured'};

% Display tolerance
data(match_rows, 9) = {expected_coords_tolerance};

% Clear any already displayed distance
data(match_rows, 10) = {[]};

% display new data with new headings
set(handles.coords_table,'Data', data);
handles.coords_table.ColumnName(5:10) = {'X expected', 'Y expected', 'Z expected', 'Is expected', 'Tolerance', 'Distance'};
handles.coords_table.ColumnWidth(5:7) = handles.coords_table.ColumnWidth(4);
handles.coords_table.ColumnWidth(8) = {'auto'};
handles.coords_table.ColumnWidth(9:10) = handles.coords_table.ColumnWidth(4);



function handles = remove_expected_coords(handles)
% handles    structure with handles and user data (see GUIDATA)
if size(handles.coords_table.Data, 2) > 1
    handles.coords_table.Data = handles.coords_table.Data(:, 1:4);
end
handles.coords_table.ColumnName = handles.coords_table.ColumnName(1:4);



% --------------------------------------------------------------------
function menu_options_expected_coords_Callback(hObject, eventdata, handles)
% hObject    handle to menu_options_expected_coords (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% disable measurements
handles.disable_measurements = true;
guidata(hObject,handles);

% Find interface objects that are set to 'on' i.e. enabled...
InterfaceObj=findobj(handles.figure1,'Enable','on');
% ... and turn them off.
set(InterfaceObj,'Enable','off');

%--------------------DIALOGUE BOX-----------------------
choice = questdlg('OPTION 1/2: Save expected coordinates between sessions?', ...
	'Expected coordinates option 1/2', ...
	'Yes','No','No');
switch choice
    case 'Yes'
        handles.save_expected_coords = true;
    if ~isfield(handles, 'expected_coords')
        uiwait(warndlg('No expected coordinates have been loaded yet. If the program is restarted without loading expected coordinates, a warning will be emitted unless "No" is selected.', 'Warning'))
    end
    case 'No'
        handles.save_expected_coords = false;
end
choice = questdlg('OPTION 2/2: Disallow measurements until expected coordinates loaded?', ...
	'Expected coordinates option 2/2', ...
	'Yes','No','No');
switch choice
    case 'Yes'
        handles.only_measure_with_expected_coords = true;
    case 'No'
        handles.only_measure_with_expected_coords = false;
end

% Re-enable the interface objects.
set(InterfaceObj,'Enable','on');

% re-enable measurements
handles.disable_measurements = false;
guidata(hObject,handles);


% --- Executes on button press in pushbutton_import_expected_coordinates.
function pushbutton_import_expected_coordinates_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_import_expected_coordinates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
menu_file_import_expected_coordinates_Callback(hObject, eventdata, handles)


% --- Executes on button press in pushbutton_insert_rows_below.
function pushbutton_insert_rows_below_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_insert_rows_below (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% disable measurements
handles.disable_measurements = true;
guidata(hObject,handles);

% Find interface objects that are set to 'on' i.e. enabled...
InterfaceObj=findobj(handles.figure1,'Enable','on');
% ... and turn them off.
set(InterfaceObj,'Enable','off');

if strcmp(handles.menu_options_insert_by_number.Checked, 'on')
    assert(strcmp(handles.menu_options_insert_by_location_name.Checked, 'off'));
    prompt = {'Number of rows:'};
    dlg_title = 'Insert Rows Below...';
    num_lines = 1;
    if isfield(handles,'previousInsertRowsVal')
        defaultans = {num2str(handles.previousInsertRowsVal)};
    else
        defaultans = {'1'};
    end
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    if isempty(answer)
        % cancel selected

        % Re-enable the interface objects.
        set(InterfaceObj,'Enable','on');

        % re-enable measurements
        handles.disable_measurements = false;
        guidata(hObject,handles);
        return
    end
    try
        nrows = str2num(answer{1});
        assert(nrows > 0);
        assert(floor(nrows) == nrows);
        handles.previousInsertRowsVal = nrows;
    catch
        h = errordlg('Number of rows to insert must be a positive integer.', 'Error', 'modal');
        uiwait(h)
        % Re-enable the interface objects.
        set(InterfaceObj,'Enable','on');
        % re-enable measurements
        handles.disable_measurements = false;
        guidata(hObject,handles);
        return
    end
else
    assert(strcmp(handles.menu_options_insert_by_location_name.Checked, 'on'))
    assert(strcmp(handles.menu_options_insert_by_number.Checked, 'off'))
    % Open pre-populate rows dialogue
    valid_entries = false;
    while ~valid_entries
        prompt = {'Start number (leave blank to not use):', 'End number (leave blank to not use):', 'Start letter (leave blank to not use):', 'End letter (leave blank to not use):'};
        dlg_title = 'Names to use (e.g. 1a, 1b, 2a, 2b...)';
        num_lines = 1;
        defaultans = {'', '', '', ''};
        if isfield(handles,'previousStartNumber')
            defaultans{1} = num2str(handles.previousStartNumber);
        end
        if isfield(handles,'previousEndNumber')
            defaultans{2} = num2str(handles.previousEndNumber);
        end
        if isfield(handles,'previousStartLetter')
            defaultans{3} = num2str(handles.previousStartLetter);
        end
        if isfield(handles,'previousEndLetter')
            defaultans{4} = num2str(handles.previousEndLetter);
        end
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        if isempty(answer)
            % cancel selected - don't save handles, just return

            % Re-enable the interface objects.
            set(InterfaceObj,'Enable','on');

            % re-enable measurements
            handles.disable_measurements = false;

            return
        end
        try
            start_number = str2num(answer{1});
            if ~isempty(answer{1})
                assert(~isempty(start_number));
            end
            if ~isempty(start_number)
                assert(floor(start_number) == start_number);
            end
            handles.previousStartNumber = start_number;
        catch
            h = errordlg('Start number must be an integer!', 'Error', 'modal');
            uiwait(h)
            continue
        end
        try
            end_number = str2num(answer{2});
            if ~isempty(answer{2})
                assert(~isempty(end_number));
            end
            if ~isempty(end_number)
                assert(floor(end_number) == end_number);
            end
            handles.previousEndNumber = end_number;
        catch
            h = errordlg('End number must be an integer!', 'Error', 'modal');
            uiwait(h)
            continue
        end
        if isempty(start_number) && ~isempty(end_number)
            h = errordlg('If start number is unspecified, end number must also be unspecified.', 'Error', 'modal');
            uiwait(h)
            continue
        end
        if ~isempty(start_number) && isempty(end_number)
            h = errordlg('If start number is specified, end number must also be specified.', 'Error', 'modal');
            uiwait(h)
            continue
        end
        if end_number < start_number
            h = errordlg('End number must be greater than or equal to the start number!', 'Error', 'modal');
            uiwait(h)
            continue
        end
        try
            start_letter = char(answer{3});
            assert(length(start_letter) < 2);
            handles.previousStartLetter = start_letter;
        catch
            h = errordlg('Only one start letter can be specified.', 'Error', 'modal');
            uiwait(h)
            continue
        end
        try
            end_letter = char(answer{4});
            assert(length(end_letter) < 2);
            handles.previousEndLetter = end_letter;
        catch
            errordlg('Only one end letter can be specified.', 'Error', 'modal');
            continue
        end
        if isempty(start_number) && ~isempty(end_number)
            h = errordlg('If start letter is unspecified, end letter must also be unspecified.', 'Error', 'modal');
            uiwait(h)
            continue
        end
        if ~isempty(start_number) && isempty(end_number)
            h = errordlg('If start letter is specified, end letter must also be specified.', 'Error', 'modal');
            uiwait(h)
            continue
        end
        if end_letter < start_letter
            h = errordlg('End letter must be after or the same as the start letter!', 'Error', 'modal');
            uiwait(h)
            continue
        end
        valid_entries = true;
        % work out nrows to add
        if isempty(start_number) && ~isempty(start_letter)
            nrows = end_letter - start_letter;
        elseif ~isempty(start_number) && isempty(start_letter)
            nrows = end_number - start_number;
        else
            nrows = (1 + end_letter - start_letter)*(1 + end_number - start_number);
        end
    end
end

if isempty(nrows)
    h = warndlg('No rows have been specified!', 'Insert Rows Below...');
    uiwait(h)
    % Re-enable the interface objects.
    set(InterfaceObj,'Enable','on');
    % re-enable measurements
    handles.disable_measurements = false;
    guidata(hObject,handles);
    return
end

[error, handles] = insert_rows_below(hObject, eventdata, handles, nrows);

if strcmp(handles.menu_options_insert_by_location_name.Checked, 'on') && ~error
    assert(strcmp(handles.menu_options_insert_by_number.Checked, 'off'))
    % Fill in names
    assert(isfield(handles,'selectedRow'))
    data = get(handles.coords_table,'Data');
    start_row = handles.selectedRow(end) + 1;
    end_row = start_row + nrows - 1;
    current_number = start_number;
    current_letter = start_letter;
    for i = start_row:end_row
        % loop through letters, then through numbers
        name = strcat(num2str(current_number), current_letter);
        reset_letter = current_letter >= end_letter;
        if ~isempty(reset_letter) && reset_letter
            current_letter = start_letter;
            reset_number = current_number >= end_number;
            if ~isempty(reset_number) && reset_number
                current_number = start_number;
            else
                current_number = current_number + 1;
            end
        else
            current_letter = char(current_letter + 1);
        end
        data{i,1} = name;
    end
    set(handles.coords_table,'Data', data)
end

% Re-enable the interface objects.
set(InterfaceObj,'Enable','on');

% re-enable measurements
handles.disable_measurements = false;
guidata(hObject,handles);


function [error, handles] = insert_rows_below(hObject, eventdata, handles, nrows)

error = false;

data = get(handles.coords_table,'Data');

% See if the selectedRow variable exists within the handles struct
% (doesn't if no selection performed before clicking or cell has been
% deselected)
if(isfield(handles,'selectedRow'))
    if(handles.selectedRow(end) < size(handles.AtlasLandmarks, 1))
        errordlg('Cannot insert or delete Atlas Points','Error','modal');
        error = true;
        return
    else
        row = handles.selectedRow(end);
        % insert new rows
        data = [data(1:row,:); cell(nrows,size(data,2)); data(row+1:end,:)];

        % check if have added rows within where measurement has already been
        % made
        if(handles.selectedRow(end) < handles.point_count)
            % increment point count to account for extra points
            handles.point_count = handles.point_count + nrows;
        end

        % Locations list has now been edited so change bool.
        handles.editedLocationsList = true;

    end
else
%    % insert empty row at the end
%    data{end+1,1} = [];

    % Tell user to select a row before inserting
    errordlg('Please select a row to insert below.','Insert Error','modal');
    error = true;
    return
end
% save the newly changed data to the table on the gui
set(handles.coords_table,'Data',data);
guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_edit_insert_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_insert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_edit_delete_rows_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_delete_rows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DeleteRowPushbutton_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_edit_insert_rows_below_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_insert_rows_below (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pushbutton_insert_rows_below_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_edit_insert_1_row_below_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_insert_1_row_below (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
InsertRowPushbutton_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_file_reset_all_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_reset_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
choice = questdlg('Are you sure you want to reset? Any unsaved measurements will be lost! Settings will be retained.', ...
	'Are you sure?', ...
	'Yes','No','No');
switch choice
    case 'Yes'
        % clear previous measurements and headmap from plot...
        cla(handles.coord_plot);
        % and replot axes.
        axis(handles.coord_plot,'equal');
        save_settings(handles);
        save_locations(handles);
        handles = close_serial_port(handles);
        DIGIGUI_OutputFcn(hObject, eventdata, handles);
    case 'No'
        return;
end


% --------------------------------------------------------------------
function menu_file_reset_expected_coords_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_reset_expected_coords (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'expected_coords')
    handles = rmfield(handles, 'expected_coords');
end
if isfield(handles, 'expected_coords_tolerance')
    handles = rmfield(handles, 'expected_coords_tolerance');
end
handles = remove_expected_coords(handles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_options_location_prepopulation_Callback(hObject, eventdata, handles)
% hObject    handle to menu_options_location_prepopulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% disable measurements
if strcmp(handles.menu_options_location_prepopulation.Checked, 'on')
    handles.menu_options_location_prepopulation.Checked = 'off';
else
    handles.menu_options_location_prepopulation.Checked = 'on';
end
guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_options_insert_Callback(hObject, eventdata, handles)
% hObject    handle to menu_options_insert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_options_insert_by_number_Callback(hObject, eventdata, handles)
% hObject    handle to menu_options_insert_by_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.menu_options_insert_by_number.Checked, 'on')
    handles.menu_options_insert_by_number.Checked = 'off';
    handles.menu_options_insert_by_location_name.Checked = 'on';
else
    handles.menu_options_insert_by_number.Checked = 'on';
    handles.menu_options_insert_by_location_name.Checked = 'off';
end
guidata(hObject,handles);

% --------------------------------------------------------------------
function menu_options_insert_by_location_name_Callback(hObject, eventdata, handles)
% hObject    handle to menu_options_insert_by_location_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.menu_options_insert_by_location_name.Checked, 'on')
    handles.menu_options_insert_by_number.Checked = 'on';
    handles.menu_options_insert_by_location_name.Checked = 'off';
else
    handles.menu_options_insert_by_number.Checked = 'off';
    handles.menu_options_insert_by_location_name.Checked = 'on';
end
guidata(hObject,handles);
