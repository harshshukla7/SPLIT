function writeData(name, x, type)
%
% Write data to the file as a constant vector
%

persistent data
if isempty(data)
  data = struct;
end

if nargin == 0 || strcmpi(name, 'clear')
  data = struct;  
  return
end

if nargin == 1
  f = fopen('name', 'w');
  
  % Write the header out
  hdr = {sprintf(['Optimization data created by SPLIT\n'...
    '%s\n',...
    'More information on SPLIT toolbox: la.epfl.ch'], datestr(datetime)),...
    length(data)};
  
  fclose(f);
else
  i = length(data)+1;
  data(i).name = name;
  data(i).x    = x;
  data(i).type = type;
end


fwrite(dataFile, uint32(length(x(:))), 'uint32');
name = sprintf('%-20s',name);
fwrite(dataFile, name, 'char');
% fwrite(dataFile, 


% if nargin < 3, type = 'double'; end
% 
% str = sprintf(' %.20g,', x);
% str(end) = [];
% 
% p('static const %s %s[] = {%s};', type, name, str);


% sep = '';
% for i = 1:length(x)
%   p('%s%.20g', sep, x(i));
%   sep = ', ';
% end
% pl('};')
