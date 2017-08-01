function p(varargin)

global f

if length(varargin) == 0
  return
end

% Test if the first argument is a file handle
fNum = f;
if isnumeric(varargin{1}) || isscalar(varargin{1})
  a = fopen(varargin{1});
  if ~isempty(a) & (a ~= -1)
    fNum = varargin{1};
    varargin = {varargin{2:end}};
  end
end

fprintf(fNum, varargin{:});
