function varargout = BORL_poldigi_gui_ver1(varargin)
% -------------------------------------------------------------------------
% Copyright (C) Reuben William Hill 2013, working at BORL, UCL, London
% reuben.w.hill@gmail.com
%
% -------------------------------------------------------------------------
%                   -- BORL POLGUI v1.1.0 for Matlab R2015b  -- 
% -------------------------------------------------------------------------
%
% For the Polhemus PATRIOT digitiser, attached to stylus pen with button.
%
% A list of points to digitise is imported from a text file, where each 
% point is on a new line.
%
% The baud rate is set via the variable "BaudRate" in 
% BORL_poldigi_gui_ver1_OpeningFcn and has default value 115200.
%
% Points are digitised by pressing the stylus button.
%
% After getting 5 reference cardinal points 'Nasion','Inion','Ar','Al' and 
% 'Cz', an Atlas reference baby head, with these points marked, is mapped
% onto the graph display of points.
%
% Before the allignment of the head model, a coordinte transform is done: 
% 1: place the 'inion' at the origin
% 2: rotate the 'Al' into the y axis
% 3: rotate the 'Ar' into the xy plane about the new 'Inion'-'Al' y axis
% 4: Rotate about the 'Al'-'Ar' axis to bring the 'Nasion' into the xy
%    plane, thus alligning the inion and nasion.
% This coordinate transform is then applied to all measured points.
%
% A tab delimited list of points and their XYZ coords is outputted to 
% a file of the users choosing.
%
%
% MATLAB GUIDE Generated comments:
% BORL_POLDIGI_GUI_VER1 MATLAB code for BORL_poldigi_gui_ver1.fig
%      BORL_POLDIGI_GUI_VER1, by itself, creates a new BORL_POLDIGI_GUI_VER1 or raises the existing
%      singleton*.
%
%      H = BORL_POLDIGI_GUI_VER1 returns the handle to a new BORL_POLDIGI_GUI_VER1 or the handle to
%      the existing singleton*.
%
%      BORL_POLDIGI_GUI_VER1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BORL_POLDIGI_GUI_VER1.M with the given input arguments.
%
%      BORL_POLDIGI_GUI_VER1('Property','Value',...) creates a new BORL_POLDIGI_GUI_VER1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BORL_poldigi_gui_ver1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BORL_poldigi_gui_ver1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BORL_poldigi_gui_ver1

% Last Modified by GUIDE v2.5 19-May-2016 13:57:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BORL_poldigi_gui_ver1_OpeningFcn, ...
                   'gui_OutputFcn',  @BORL_poldigi_gui_ver1_OutputFcn, ...
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


% --- Executes just before BORL_poldigi_gui_ver1 is made visible.
function BORL_poldigi_gui_ver1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BORL_poldigi_gui_ver1 (see VARARGIN)

% Choose default command line output for BORL_poldigi_gui_ver1
handles.output = hObject;

%--------------------define close request function----------------------
%the function "CloseFcn" that I define now runs when quitting the gui
set(gcf,'CloseRequestFcn',{@CloseFcn,handles});

%------------------------CREATE SERIAL OBJECT---------------------------

% Create serial object and set baud rate 
BaudRate = 115200;

handles.COMport = FindPatriotSerial(BaudRate);

if(handles.COMport ~= 0) %patriot found
    handles.serial = serial(handles.COMport,'BaudRate', BaudRate);
else
%----------------------ERROR IF DEVICE NOT FOUND------------------------   
    str1 = 'Polhemus Patriot Device not found or communicated with successfully.';
    str2 = ['Check the device is on and its baud rate is set to ' sprintf('%i',BaudRate) '.'];
    str3 = 'If running in MATLAB, try restarting MATLAB to scan for new serial devices.';
    str4 = 'Also consider turning the device off and on. Take care to give the device time to reinitialise before trying again.';
    errstr = sprintf('%s\n\n%s\n\n%s\n\n%s',str1,str2,str3,str4);
    %display error message
    uiwait(errordlg(errstr,'Polhemus Communications Initialisation Error'));
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BORL_poldigi_gui_ver1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);





% --- Outputs from this function are returned to the command line.
function varargout = BORL_poldigi_gui_ver1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


if(handles.COMport == 0) %quit if com port not found
    CloseFcn(hObject,eventdata,handles);
else
    %------------------------Initialise Variables-----------------------
    %Set the initial point count to 0. This is incremented before each
    %measured head point until the last head point is measured.
    handles.point_count = 0;

    handles.all_points_found = false;

    %--------------------HEADPOINTS TO DIGITISE INPUT-----------------------

    [filename,pathname] = ... 
        uigetfile({'*.txt;*.dat;*.csv','Text Files (*.txt) (*.dat) (*.csv)'} ...
                  ,'Select Location List File - Each Measurement Point Should be on a New Line');

    if isequal(filename,0)
        disp('User selected Cancel')
        %Quit the gui
        CloseFcn(hObject,eventdata,handles);

    else % the below code runs if a file is selected...

        %load data needed for headpoint plotting
        handles.AtlasLandmarks = load('refpts_landmarks.mat');
        handles.AtlasLandmarks = handles.AtlasLandmarks.pts;
        handles.mesh = load('scalpSurfaceMesh.mat');
        handles.mesh = handles.mesh.mesh;

        disp(['User selected ', fullfile(pathname, filename)])

        % read headpoints text file
        
%        if isdeployed
%           addpath pathname;
%        end
        
        FileID = fopen([pathname filename]);
        locations = textscan(FileID,'%s','delimiter','\n');
        % locations is a local variable that holds location data in this
        % function

        % append to list of reference points and convert to string array
        locations = ['Nasion';'Inion';'Ar';'Al';'Cz'; ... 
                                                locations{1,1}];
        fclose(FileID);

        %error test the first serial port functions...
        try 
            %------------------------SERIAL CALLBACK SETUP---------------------
            %setup callback function to run when the polhemus system sends 48 bytes
            %48 bytes is generally position data
            handles.serial.BytesAvailableFcnCount = 48;
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
            disp('COM PORT ERROR OCCURRED: Check COM1 Connection. Baud rate should be 115200')
            %run close function to close gui and delete serial port objects if
            %error occurs.
            CloseFcn(hObject,eventdata,handles);
            error(message('MATLAB:serial:fopen:opfailed', serialException.message))
        end

        %-------------set graph axes labels and properties---------------------

        xlabel(handles.coord_plot,'X');
        ylabel(handles.coord_plot,'Y');
        zlabel(handles.coord_plot,'Z');
        
        guidata(hObject, handles);
        
    end
   
end
% Update handles structure



% --- Executes on button press in HeadAlign.
function HeadAlign_Callback(hObject, eventdata, handles)
% hObject    handle to HeadAlign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.point_count >= 5)
    
    % extract the locations
    locations = handles.coords_table.Data;
    
    % extract the landmark locations (the first five data points)...
    landmarks = locations(1:5,2:4);
    % ... and convert to ordinary array from cell array
    landmarks = cell2mat(landmarks);
    
    %get transformation matrix to new coord system
    [TransformMatrix,TransformVector] = GetCoordTransform(landmarks);
    %save tranformation
    handles.TransformMatrix = TransformMatrix;
    handles.TransformVector = TransformVector;
    
    % reset list of points to just show locations to find so transformed 
    % points can be plotted
    locations = handles.coords_table.Data;
    
    hold on
    
    for k = 1:size(landmarks,1)
        %transform cardinal points
        landmarks(k,:) = landmarks(k,:) + TransformVector;
        landmarks(k,:) = landmarks(k,:)*TransformMatrix';

        %remove old point from graph
        delete(handles.pointhandle(k));
        
        %replot point
        handles.pointhandle(k) = plot3(landmarks(k,1), ...
                                       landmarks(k,2), ...
                                       landmarks(k,3), ...
                                       'm.', 'MarkerSize', 20, ...
                                       'Parent' , handles.coord_plot);
        
        %replot axes...
        axis(handles.coord_plot,'equal');
        
        %update newly transformed cardinal point coords (converting back
        %to a cell array first)
        locations(k,2:4) = num2cell(landmarks(k,1:3));
                                              
    end
    
    hold off
    
    
    % Show newly transformed cardinal point coords on table
    handles.coords_table.Data = locations;
    
    %find matrix (A) and vector (B) needed to map head to cardinal points
    %with affine transformation
    [A,B] = affinemap(handles.AtlasLandmarks,landmarks);

    mesh_trans = handles.mesh;
    mesh_trans.node = affine_trans_RJC(handles.mesh.node,A,B);

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
    light;
    lighting gouraud;
    axis equal;
    hold off;
    
    
    % disable  headalign button
    set(hObject,'Enable','off');
    
    % Update handles structure
    guidata(hObject, handles);

end


%This function outputs a vector to transform inion to origin
%also outputs rotation transformation matrix to allign the head model to 
%intuitive coordinates. 
%To apply to row vector, the vector should be multiplied by the transpose
%with the vector on the left
function [Matrix,vector] = GetCoordTransform(landmarks)
%calculate lengths between vectors

%these are the untransformed reference points, defined here in case I want
%to use them
Nasion = landmarks(1,:);
Inion = landmarks(2,:);
Ar = landmarks(3,:);
Al = landmarks(4,:);
Cz = landmarks(5,:);

%------TRANSLATION-------

%translate inion to origion
vector = -Inion;

%translate Al
Al = Al + vector;

%------ROTATE AL TO Y AXIS-------

%calculate rotation to y axis
AlToYAxisRot = vrrotvec(Al,[0,1,0]);

%convert to rotation matrix
AlToYAxisMatrix = vrrotvec2mat(AlToYAxisRot);

%repmat
% Apply translation and rotation to points
for k = 1:5
    landmarks(k,:) = landmarks(k,:)+vector;
    landmarks(k,:) = landmarks(k,:)*AlToYAxisMatrix';
end

%------ROTATE AR TO XY PLANE ABOUT INION-AL AXIS-------

% find angle to rotate nasion into XY plane about the new y axis
[ArToXYRotAngle,~] = cart2pol(landmarks(3,1),landmarks(3,3));

%Find second rotation matrix
ArToXYMatrix = vrrotvec2mat([0,1,0,ArToXYRotAngle]);

% Apply second rotation to points
for k = 1:5
    landmarks(k,:) = landmarks(k,:)*ArToXYMatrix';
end

%------FINAL ROTATION ABOUT AL-AR AXIS TO ALLIGN INION AND NASION-------

%find angle of nasion to xy plane
[~,NasionToXYRotAngle,~] = cart2sph(landmarks(1,1),landmarks(1,2),landmarks(1,3));
%define vector to rotate around (the line joining AL and AR)
NasionRotVector = landmarks(4,:) - landmarks(3,:);
%find rotation matrix
NasionToXYRotMatrix = vrrotvec2mat([NasionRotVector, NasionToXYRotAngle]);

%------OUTPUT FINAL MATRIX-------
Matrix = NasionToXYRotMatrix*ArToXYMatrix*AlToYAxisMatrix;





function CloseFcn(source,event,handles)
%my user-defined close request function
%closes the serial port

handles = guidata(handles.figure1);

%close port only if not closed
if(handles.COMport ~= 0)
    if( ~ strcmp(handles.serial.status, 'closed') )
        fclose(handles.serial);
    end
end

delete(gcf);


function ReadCoordsCallback(s,BytesAvailable,handles)

% Update handles structure to most current version
handles = guidata(handles.figure1);

%don't run most of the callback if all points have been found..
if(~handles.all_points_found)
    
    %don't run most of the callback if waiting to do alignment...
    if(handles.point_count < 5 || strcmp(get(handles.HeadAlign,'Enable'),'off') )

        %increment the point count before measurement
        handles.point_count = handles.point_count + 1;

        %read the data on the serial port that triggered the callback
        data_str=fgetl(s);

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
        Coords = data_num(1,2:4);

        % disable head alignment butten for first five points (they are the
        % landmark positions)
        if(handles.point_count < 5)
            set(handles.HeadAlign,'Enable','off');
        % enable head allign after 5 points...
        elseif(handles.point_count == 5)
            set(handles.HeadAlign,'Enable','on');
        % Do coord transform on points measured after landmark points
        else 
            Coords = Coords + handles.TransformVector;
            Coords = Coords*handles.TransformMatrix';
        end

        % Update table with newly measured x y and z values
        handles.coords_table.Data(handles.point_count,2:4) = ...
            num2cell(Coords);
        
        % update point to look for (unless at end of list as given by the
        % length of handles.coords_table.Data - ie the number of headpoints)
        if( handles.point_count < size(handles.coords_table.Data,1) )
            set(handles.infobox,'string',...
                handles.coords_table.Data(handles.point_count+1,1));
                % (Set to the next position on the table)
        else
            set(handles.infobox,'string','All points Collected!');
            handles.all_points_found = true;
        end  

        %add the measured point to the 3d graph
        hold(handles.coord_plot,'on');
        %save the handle of the point so it can be removed later...
        if(handles.point_count <= 5)
            handles.pointhandle(handles.point_count) = plot3(Coords(1), ...
                                                Coords(2),Coords(3), ...
                                                'm.', 'MarkerSize', 20, ...
                                                'Parent' , handles.coord_plot);
        else %Note: above marker points are plotted differently
            handles.pointhandle(handles.point_count) = plot3(Coords(1), ...
                                                Coords(2),Coords(3), ... 
                                               'b.', 'MarkerSize', 20, ...
                                               'Parent',handles.coord_plot);
        end
        hold(handles.coord_plot,'off'); 
        %replot axes...
        axis(handles.coord_plot,'equal');
    end
end

% Update handles structure
guidata(handles.figure1,handles);


% --- Executes on button press in remove_last_pt.
function remove_last_pt_Callback(hObject, eventdata, handles)
% hObject    handle to remove_last_pt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%don't delete points if alignment already done or at first point
if (handles.point_count ~= 0)
    if(handles.point_count ~= 5 || strcmp(get(handles.HeadAlign,'Enable'),'on') )

        % Set the last measured values of x, y and z to be empty cells
        handles.coords_table.Data{handles.point_count,2} = []; % x
        handles.coords_table.Data{handles.point_count,3} = []; % y
        handles.coords_table.Data{handles.point_count,4} = []; % z

        % Remove point from graph...
        delete(handles.pointhandle(handles.point_count));
        % and replot axes.
        axis(handles.coord_plot,'equal');
        
        % Decrement point_count so next measurement is of the point which
        % has just been deleted
        handles.point_count = handles.point_count - 1;
        
        % Update the all points collected bool if set to true
        if(handles.all_points_found)
            handles.all_points_found = false;
        end

        % Disable align if now not enough points
        if(handles.point_count <= 5)
            set(handles.HeadAlign,'Enable','off');
        end

        % Update handles structure
        guidata(handles.figure1,handles);

    end
end



% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Open a "Save As..." Dialogue with different saving options as shown.
% The filterIndex gives the index (1, 2 or 3) of the chosen save type.
[fileName,pathName,filterIndex] = ... 
    uiputfile({'*.csv;*.dat;*.txt', ... 
    'Comma-delimited text files (*.csv) (*.dat) (*.txt)'; ...
    ...
    '*.mat','MAT-file (*.mat)'; ...
    ...
    '*.xls;*.xlsb;*.xlsm;*.xlsx', ...
    'Excel® spreadsheet files (*.xls) (*.xlsb) (*.xlsm) (*.xlsx)'; ...
    },'Save As...');

% If the chosen save type is .mat then use a standard matlab save command
if(filterIndex == 2)
    dataOutput = get(handles.coords_table,'Data');
    save([pathName fileName],'dataOutput');
% Otherwise create a table from the cell array and output that to file.
else
    tableToOutput = cell2table(get(handles.coords_table,'Data'),...
                    'VariableNames',get(handles.coords_table,'ColumnName'));
    % Note that writetable changes its output depending on the fileName
    % type.
    writetable(tableToOutput,[pathName fileName]);
end


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
    % extract the row from where the user clicked on the table.
    handles.selectedRow = eventdata.Indices(1);
end
guidata(hObject,handles)



% --- Executes on button press in InsertRowPushbutton.
function InsertRowPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to InsertRowPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = get(handles.coords_table,'Data');

% See if the selectedRow variable exists within the handles struct 
% (doesn't if no selection performed before clicking or cell has been 
% deselected)
if(isfield(handles,'selectedRow'))
    
    if(handles.selectedRow < 5)
        % trying to insert into one of the landmark points
        warningString = {'About to insert or delete a landmark point.'; ... 
                    'Pressing OK will remove the ability to align atlas points...'}
        button = questdlg(warningString,'Warning','OK','Cancel','default');
        if(button == 'OK')
            handles.CannotAlignAtlasPoints = true;
        else
            % user selected "Cancel"
            set
        end      
    end
    
    % insert below selected row...
    row = handles.selectedRow;
    dataBelowSelectedRow = data(row+1:end,:);
    % add new row by adding a single rowed cell array
    data(row+1,:) = cell(1,size(data,2));
    % add back the data that was saved before by concatenating below where
    % the new row has been added.
    data = [data(1:row+1,:) ; dataBelowSelectedRow];  
    
else
    % insert empty row at the end
    data{end+1,1} = [];
end
% save the newly changed data to the table on the gui
set(handles.coords_table,'Data',data); 

% --- Executes on button press in DeleteRowPushbutton.
function DeleteRowPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteRowPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

