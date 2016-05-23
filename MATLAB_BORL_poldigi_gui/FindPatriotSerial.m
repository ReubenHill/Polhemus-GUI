function [ COMPort ] = FindPatriotSerial(BaudRate)
%Finds the COM port of the Polhemus Patriot device when connected by RS232
%   The function returns the COM port if successful, otherwise it returns
%   0;

disp('Now looking for Polhemus Patriot device...')
disp('Existing serial port objects in MATLAB will be deleted')

%clear any existing serial ports
delete(instrfindall('Type','serial'));

%get all serial ports detected by matlab upon it's opening...
serialInfo = instrhwinfo('serial');

%loop through all found serial ports and test comms
for k = 1:size(serialInfo.AvailableSerialPorts,1)
    try
        disp(['Now testing port ' serialInfo.AvailableSerialPorts{k,1} ])
        s = serial(serialInfo.AvailableSerialPorts{k,1},'BaudRate',BaudRate);
        fopen(s);
        %write 'p', which causes stylus coordinates to be returned
        fwrite(s,'p');
        pause(0.1);
        %read coordinates...
        A = fgetl (s);
        fclose(s);
        delete(s);
        if( size(A,2) == 59 ) %a position reading is a 59 character string
            %com port of patriot successfully found!    
            COMPort = serialInfo.AvailableSerialPorts{k,1};
            disp( ['Polhemus Patriot device successfully found on ' serialInfo.AvailableSerialPorts{k,1}] )
            disp('Device returned the following string upon "p" command')
            disp(A);
            return; %return to calling function
        end
    catch serialException %if error message, the below runs
        fclose(s);
        delete(s);
    end
end

disp('Polhemus Patriot not found!')

%if this part of function is reached without returning, the patriot device
%was not communicated with successfully
COMPort = 0;

return
    

