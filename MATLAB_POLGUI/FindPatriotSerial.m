function [ COMPort, sensors ] = FindPatriotSerial(BaudRate)
%Finds the COM port of the Polhemus Patriot device when connected by RS232
%   The function returns the COM port if successful, otherwise it returns
%   0;

disp('Now looking for Polhemus Patriot device...')

%clear any existing serial ports
disp('Any existing serial port objects in MATLAB being deleted...')
delete(instrfindall('Type','serial'));

%get all serial ports detected by matlab upon it's opening...
serialInfo = instrhwinfo('serial');

%loop through all found serial ports and test comms
for k = 1:size(serialInfo.AvailableSerialPorts,1)
    try
        disp(['Now testing port ' serialInfo.AvailableSerialPorts{k,1} '...' ])
        s = serial(serialInfo.AvailableSerialPorts{k,1},'BaudRate',BaudRate);
        fopen(s);
        %write 'p', which causes stylus coordinates to be returned
        fwrite(s,'p');
        pause(0.1);
        %read coordinates...
        A = fread (s);
        fclose(s);
        delete(s);
        if( length(A) >= 59 ) %A position reading will be a string of at least 59 characters
            %com port of patriot successfully found!    
            COMPort = serialInfo.AvailableSerialPorts{k,1};
            disp( ['Polhemus Patriot device successfully found on ' serialInfo.AvailableSerialPorts{k,1}] )
            disp('Device returned the following string upon "p" command:')
            disp(char(A'));
            %determine if output was for a 1 or 2 sensor setup
            if( length(A) < 120 ) % single source (might include "patriot ready!")
                sensors = 1;
            else
                sensors = 2;
            end
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
    

