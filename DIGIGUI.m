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
% Before the allignment of the head model, a coordinte transform is done:
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
% along with a given tolerance in a space delimited text file. Click
% 'File/Import Expected Coordinates...' for more detail.
%
% A tab delimited list of points and their XYZ coords is outputted to
% a file of the user's choosing.
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

% Last Modified by GUIDE v2.5 23-Jun-2023 10:06:17

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
    load(settings_loc,'error_distance','atlas_dir','save_expected_coords');
    disp(['error_distance: ', num2str(error_distance)]);
    disp(['atlas_dir: ', atlas_dir]);
    disp(['save_expected_coords: ', num2str(save_expected_coords)]);
catch
    uiwait(errordlg('Settings missing or not found, falling back on default settings.','Settings Error'));
    try
        load('default_settings.mat','error_distance','atlas_dir','save_expected_coords');
        disp(['error_distance: ', num2str(error_distance)]);
        if ~isdeployed
            atlas_dir = fullfile(pwd, atlas_dir);
        else
            atlas_dir = fullfile(ctfroot, atlas_dir);
        end
        disp(['atlas_dir: ', atlas_dir]);
        disp(['save_expected_coords: ', num2str(save_expected_coords)]);
    catch
        uiwait(errordlg('Default settings missing or not found! Quitting.','Settings Error'));
        %Quit the gui
        guidata(hObject, handles);
        CloseFcn(hObject,eventdata,handles);
        return
    end
end

handles.error_distance = error_distance;
handles.save_expected_coords = save_expected_coords;

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
        handles.expected_coords = expected_coords;
        handles.expected_coords_tolerance = expected_coords_tolerance;
    catch
        uiwait(warndlg('Could load the expected coordinates and tolerance from saved_expected_coords.mat.',...
            'Expected Coordinates Warning','modal'));
        if isfield(handles, 'expected_coords')
            handles = rmfield(handles, 'expected_coords');
        end
        if isfield(handles, 'expected_coords_tolerance')
            handles = rmfield(handles, 'expected_coords_tolerance');
        end
    end
else
    % remove the saved file to avoid confusion
    if exist(saved_expected_coords_loc, 'file')==2
        delete(saved_expected_coords_loc);
    end
end

% Get the checkmark icon for use when checking if coordinates match
handles.checkmark_icon = imread('checkmark.tif');

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





function CloseFcn(source,event,handles)
%my user-defined close request function
%closes the serial port

try
    handles = guidata(handles.figure1);
catch
    handles = struct([]);
end

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
if isfield(handles, 'save_expected_coords')
    disp(['Saving save_expected_coords to ', settings_loc]);
    save_expected_coords = handles.save_expected_coords;
    save(settings_loc, 'save_expected_coords', '-append');
end

%close port only if not closed
if(isfield(handles,'COMport'))
    if(handles.COMport ~= 0)
        if( ~ strcmp(handles.serial.status, 'closed') )
            fclose(handles.serial);
        end
        delete(handles.serial);
    end
end

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
% Don't measure if we have a double tap warning open or unexpected
% measurement error open
if isfield(handles, 'doubleTapWarnFigure') && isvalid(handles.doubleTapWarnFigure)
    return
end
if isfield(handles, 'unexpectedMeasurementErrorFigure') && isvalid(handles.unexpectedMeasurementErrorFigure)
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

% Check if table is currently full - if it is then adding a new point
% will expand the table...
if(handles.point_count > size(data,1))
    % ... so update the bool that tracks if location names have been
    % edited. When the user saves their data they will therefore be
    % prompted to save the locations list too.
    handles.editedLocationsList = true;
end

% Update table with newly measured x y and z values
data(handles.point_count,2:4) = num2cell(Coords);
set(handles.coords_table,'Data',data);

% Double tap warning...
if handles.point_count > 1
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
            msg = sprintf('%s measurement was only %0.2g cm from %s measurement!\nCurrent error distance is %0.2g cm. This can be changed in the options.', this_point{1}, distance, last_point{1}, handles.error_distance);
            handles.doubleTapErrorFigure = errordlg(msg, 'Double tap error', 'modal');
            % reset point count, remove data, update guidata and exit
            data(handles.point_count,2:4) = {[], [], []};
            set(handles.coords_table,'Data',data);
            handles.point_count = previous_point_count;
            guidata(handles.figure1,handles);
            return
        end
    end
end

% Remove any open expected measurement figures
if isfield(handles, 'expectedMeasurementFigure') && isvalid(handles.expectedMeasurementFigure)
    close(handles.expectedMeasurementFigure)
end
% Check against expected coordinates if any are specified
if isfield(handles, 'expected_coords')
    if handles.point_count <= size(handles.expected_coords, 1)
        distance = norm(Coords - handles.expected_coords(handles.point_count, :));
        this_point = data(handles.point_count, 1);
        if distance > handles.expected_coords_tolerance
            msg = sprintf('%s measurement is %0.2g cm from the expected location of (%0.2g, %0.2g, %0.2g) cm!\nCurrent tolerance is %0.2g cm. Tolerance is set in the expected coordinates file.', this_point{1}, distance, handles.expected_coords(handles.point_count, 1), handles.expected_coords(handles.point_count, 2), handles.expected_coords(handles.point_count, 3), handles.expected_coords_tolerance);
            handles.unexpectedMeasurementErrorFigure = errordlg(msg, 'Unexpected Measurement!', 'modal');
            % reset point count, remove data, update guidata and exit
            data(handles.point_count,2:4) = {[], [], []};
            set(handles.coords_table,'Data',data);
            handles.point_count = previous_point_count;
            guidata(handles.figure1,handles);
            return
        else
            msg = sprintf('%s measurement is within tolerance of (%0.2g, %0.2g, %0.2g) cm.\nCurrent tolerance is %0.2g cm. Tolerance is set in the expected coordinates file.', this_point{1}, handles.expected_coords(handles.point_count, 1), handles.expected_coords(handles.point_count, 2), handles.expected_coords(handles.point_count, 3), handles.expected_coords_tolerance);
            handles.expectedMeasurementFigure = msgbox(msg, 'Expected Measurement Success!', 'custom', handles.checkmark_icon);
        end
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

    % Set the last measured values of x, y and z to be empty cells
    data{last_point_count,2} = []; % x
    data{last_point_count,3} = []; % y
    data{last_point_count,4} = []; % z

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
    if(handles.point_count <= size(handles.AtlasLandmarks, 1))
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
        save([pathName fileName],'dataOutput');
        disp('Data is stored in cell array "dataOutput"');

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

            tableToOutput = cell2table(data,'VariableNames', ...
                                       get(handles.coords_table,'ColumnName'));
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

        % check if have deleted any rows where measurements have already
        % been made
        if(any(handles.selectedRow <= handles.point_count))

            % Remove point from graph...
            delete(handles.pointhandle(handles.point_count));
            % and replot axes.
            axis(handles.coord_plot,'equal');

            % find out how many of the selected rows are less than the
            % current point_count
            numToDecrement = nnz(handles.selectedRow <= handles.point_count);

            % decrement point count to account for number of fewer points
            handles.point_count = handles.point_count - numToDecrement;

        end
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
if ~isdeployed
    saved_location_names_loc = fullfile(pwd, 'savedLocationNames.mat');
else
    saved_location_names_loc = fullfile(ctfroot, 'DIGIGUI', 'savedLocationNames.mat');
end
disp(['Saving locations to: ', saved_location_names_loc]);
save(saved_location_names_loc, 'locations');

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
            end
        end
        % delete data on table
        data(handles.selectedRow,2:4) = cell(length(handles.selectedRow), 3);
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

% Save locations variable to be loaded next time
if ~isdeployed
    saved_location_names_loc = fullfile(pwd,'savedLocationNames.mat');
else
    saved_location_names_loc = fullfile(ctfroot, 'DIGIGUI', 'savedLocationNames.mat');
end
disp(['Saving locations to: ', saved_location_names_loc]);
save(saved_location_names_loc,'locations');

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
handles.error_distance = doubletapdialog(handles.error_distance);

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
        ,['Select expected coordinates list file - each measurement point should be'...
        ' on a new line. Tolerance should have its own line followed by two NaNs.']);
else
    [filename,pathname] = ...
        uigetfile({'*.txt;*.dat;*.csv', ...
        'Text Files (*.txt) (*.dat) (*.csv)'} ...
        ,['Select expected coordinates list file - each measurement point should be'...
        ' on a new line. Tolerance should have its own line followed by two NaNs.'],handles.userDir);
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

expected_coords = load([pathname filename]);

% Find tolerance - it should be folowed by two NaNs
[row, col] = find(isnan(expected_coords), 1, 'first');
if ~length(row) || ~length(col)
    errordlg('No NaNs found in expected coordinates file! Tolerance was therefore not found!')
    return
end
expected_coords_tolerance = expected_coords(row, col-1);
% remove the nan row
expected_coords(row, :) = [];
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
choice = questdlg('Save expected coordinates between sessions?', ...
	'Expected coordinates', ...
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
    return
end
try
    nrows = str2num(answer{1});
    assert(nrows > 0);
    assert(floor(nrows) == nrows);
    handles.previousInsertRowsVal = nrows;
catch
    errordlg('Number of rows to insert must be a positive integer.', 'Error', 'modal');
    return
end
insert_rows_below(hObject, eventdata, handles, nrows);


function insert_rows_below(hObject, eventdata, handles, nrows)

data = get(handles.coords_table,'Data');

% See if the selectedRow variable exists within the handles struct
% (doesn't if no selection performed before clicking or cell has been
% deselected)
if(isfield(handles,'selectedRow'))
    if(handles.selectedRow(end) < size(handles.AtlasLandmarks, 1))
        errordlg('Cannot insert or delete Atlas Points','Error','modal');
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
end
% save the newly changed data to the table on the gui
set(handles.coords_table,'Data',data);
guidata(hObject,handles);
