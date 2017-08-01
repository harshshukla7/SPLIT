function pl(varargin)
% Print the string to the file, with a newline


wrote = false;
if length(varargin) > 0
  if isnumeric(varargin{1}) || isscalar(varargin{1})
    a = fopen(varargin{1});
    if ~isempty(a) & (a ~= -1)
      if length(varargin) > 1
        p(varargin{:});
      end
      p(varargin{1}, '\n');
      wrote = true;
    end
  end
end

if ~wrote
  p(varargin{:});
  p('\n');
end
