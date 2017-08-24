function count = struct_write(~, fid,format,data)
%% A function to write a C structure to a binary file.
% struct_write simplifies the writing of a C structure to a binary file by
% providing a functionality similar to the Python 'struct' module. The
% element types of the structure are specified as a series of characters in
% a format string. struct_write then converts values passed in via a cell
% array to these binary formats and writes the them to a file handle.
%
% ARGUMENTS:
%
%   fid         File handle of an file open for writing.
%   format      Format string specifying the C structure. (see below)
%   data        cell array whos values are written to the file
%
% FORMAT SPECIFIERS:
% 
% Specifier     Bytes:        Notes
% c             2             char returned as a string character
% C             2             char returned as a number (double)
% b             2             signed char
% B             2             unsigned char
% h             2             short
% H             4             unsigned short
% i             4             int 
% I             4             unsigned int
% l             4             long
% L             4             unsigned long
% q             8             long long
% f             4             float
% d             8             double
% s             char[]        char[] returns character string
% r                           real 
%
% Note that unlike struct_read, multipliers (e.g 8c to read 8 chars) are
% not specified for struct_write. Struct_write instead takes advantage of
% the MATLAB's fwrite function which writes each element in the passed
% vector with the specified precision. Therefore one need only have a
% one-to-one correspondence between each cell entry in 'data' and each
% character specifier in the format string. Lets look at an example:
%
% EXAMPLE:
%
% Given the following C structure definition for a depth measurement from an
% echo-sounder:
%
% typedef struct 
% { 
%     float     depth1; //Depth1 
%     char      units;  // m for meters, f for feet
%     double    time_stamp; //Time stamp, number of seconds elapsed since midnight (00:00:00), January 1, 1970 
%     char[256] navstring; // GPS Position String
% } ECHO; 
%
% First create a cell array with the desired values for each variable
% in the struct:
%
% D:
% { [23.456],'m',[1242238192.870432012], ...
%  '$GPGGA,182503.00,7121.11801,N,15651.53926,W,1,12,0.8,20.49,M,-048,M,,*5
%  7' }
%
% Then open the file and write the information like this:
% fid = fopen('mybinaryfile','w');
% D = struct_write(fid,'fcd256c',D)
%
%
% Val Schmidt
% University of New Hampshire
% Center for Coastal and Ocean Mapping
% 2009
%
%%
% Questions?
% Should we return the number of bytes written rather than the number of
% items? 
% 

if length(format) ~= length(data)
    disp('ERROR, Format string is not equal to the number of data types.')
    return
end
types = format;

count=0;
for i = 1:length(types)
    
    switch types(i)
        
        case 'c'
            % char
           cnt = fwrite(fid,data{i},'char');
            
        case 'b'
            % signed char
            cnt = fwrite(fid,data{i},'schar');
            
        case 'B'
            % unsigned char
            cnt = fwrite(fid,data{i},'uchar');
        case 'h'
            % short
            cnt = fwrite(fid,data{i},'short');
            
        case 'H'
            % unsigned short
            cnt = fwrite(fid,data{i},'ushort');
            
        case 'i'
            % int
            cnt = fwrite(fid,data{i},'int');
            
        case 'I'
            % unsigned int
            cnt = fwrite(fid,data{i},'uint');
            
        case 'l'
            % long
            cnt = fwrite(fid,data{i},'long');
            
        case 'L'
            % unsigned long
            cnt = fwrite(fid,data{i},'ulong');
            
        case 'q'
            % long long
            
        case 'Q'
            % unsigned long long
            
        case 'f'
            % float
            cnt = fwrite(fid,data{i},'float');
        case 'd'
            % double
            cnt = fwrite(fid,data{i},'double');
            
        case 's'
            % char[]  ??? does this require adding 1 to mult for \null?
            cnt = fwrite(fid,data{i},'char');
        
        case 'r'
            % real, only for split implementation
            cnt = fwrite(fid,data{i},'char');
    end
    count = count + cnt;
end

    

end
