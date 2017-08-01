function str = genSparseMatVec(cdr,A,x,y,Ix,Iy,returnStr,tol)
%
% str = genSparseProduct(A,x,y,Ix,Iy)
%
% Generate matlab code to multiply y(Iy) = A*x(Ix), where A is very sparse
%
% A = constant matrix
% x = string of the 'x' variable name
% y = string of the 'y' variable name (if empty and A is a vector, then
% doesn't generate y(i) = part or a semi-colon at the end)
% Ix,Iy = index sets (if empty they default to [1:n])
%
% Returns cell-array of strings, unless there's only one row - then it's
% just a string
%

INNER_PRODUCT_UNROLL_LIMIT = 10;

if nargin < 8, tol = 1e-5; end
if nargin < 7, returnStr = false; end
if nargin < 6, Iy = []; end
if nargin < 5, Ix = []; end
if isempty(Iy), Iy = 1:size(A,1); end
if isempty(Ix), Ix = 1:size(A,2); end
if nargin < 7, returnStr = false; end
if nargin < 4, y = 'y'; end

if size(A,1) > 1 && size(A,2) > 1
  if isempty(y)
    y = 'y';
    warning('Defaulting to y = ''y'' since A is not a vector')
  end
end
genY = ~isempty(y);


str = {};
% if genY
%   str{end+1} = sprintf('y = coder.nullcopy(zeros(%i,1));', size(A,1));
% end

for i = 1:size(A,1)
  I = find(abs(A(i,:)) > tol);
  if isempty(I)
    str{end+1} = '+0.0';
  else
    if length(I) == 1
      a = A(i,I);
      if a == 1
        str{end+1} = sprintf('+%s(%i)', x, Ix(I));
      elseif a == -1
        str{end+1} = sprintf('-%s(%i)', x, Ix(I));
      else
        str{end+1} = sprintf('+%f*%s(%i)', a, x, Ix(I));
      end
    else
      if length(I) < INNER_PRODUCT_UNROLL_LIMIT
        t = [];
        for j = 1:length(I)
          a = A(i,I(j));
          if a == 1
            t = sprintf('%s + %s(%i)', t, x, Ix(I(j)));
          elseif a == -1
            t = sprintf('%s - %s(%i)', t, x, Ix(I(j)));
          else
            t = sprintf('%s + %f*%s(%i)', t, a, x, Ix(I(j)));
          end
        end
        str{end+1} = t;
      else
        Istr = sprintf('%i,', Ix(I));
        a = sprintf('%f,', A(i,I));
        str{end+1} = sprintf('+[%s]*%s([%s])', a(1:end-1), x, Istr(1:end-1));
      end
    end
  end
  if genY, str{end} = sprintf('%s(%i) = %s;', y, Iy(i), str{end}); end
end

if nargout == 0
  for i = 1:length(str)
    cdr.print([str{i} '\n']);
  end
else
  if returnStr
    if length(str) == 1
      str = str{1};
    end
    if length(str) == 0
      str = '';
    end
  end
end
