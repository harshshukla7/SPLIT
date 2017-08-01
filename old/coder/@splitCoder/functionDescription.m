function functionDescription(cdr, f, prefix)
% Overload display
%
% f      = file handle to print the description to. default = 1
% prefix = string to add to the front of every line. default = ''

if nargin < 3, prefix = ''; end
if nargin < 2, f = 1; end

fprintf(f, '%s', prefix); fprintf(f, '\n');
fprintf(f, '%s', prefix); fprintf(f, ' Template file:             %s\n',     cdr.templateName);
fprintf(f, '%s', prefix); fprintf(f, ' Generated matlab function: %s\n',     cdr.fname);
fprintf(f, '%s', prefix); fprintf(f, ' Generated mex function:    %s_mex\n', cdr.fname);

fprintf(f, '%s', prefix); fprintf(f, '\n');
fprintf(f, '%s', prefix); fprintf(f, ' Problem:\n');
fprintf(f, '%s', prefix); fprintf(f, '  %i variables\n', size(cdr.A,2));
fprintf(f, '%s', prefix); fprintf(f, '  %i parameters\n', size(cdr.pF,2));
fprintf(f, '%s', prefix); fprintf(f, '  %i proximal functions:\n', length(cdr.proxInd)-1);
for i = 1:length(cdr.prob.prox)
  fprintf(f, '%s', prefix); fprintf(f, '    - %s\n', cdr.prob.prox(i).typeName);
end

fprintf(f, '%s', prefix); fprintf(f, '\n');
fprintf(f, '%s', prefix); fprintf(f, ' Function prototype:\n');
fprintf(f, '%s', prefix); fprintf(f, '  [%s] = %s(%s)\n', cdr.csvList({cdr.outArgs.name}), cdr.fname, cdr.csvList({cdr.inArgs.name}));
fprintf(f, '%s', prefix); fprintf(f, '  [%s] = %s_mex(%s)\n', cdr.csvList({cdr.outArgs.name}), cdr.fname, cdr.csvList({cdr.inArgs.name}));
fprintf(f, '%s', prefix); fprintf(f, '\n');

dat(1).desc  = 'Arguments to the generated function';
dat(1).dat   = cdr.inArgs;
dat(1).defaultNum = 0;
dat(2).desc  = 'Outputs of the generated function';
dat(2).dat   = cdr.outArgs;
dat(2).defaultNum = 0;
dat(3).desc  = 'Parameters available at compile-time and run-time';
dat(3).dat   = cdr.inConstArgs;
dat(3).defaultNum = 11;
dat(4).desc  = 'Parameters available at compile-time, but not run-time';
dat(4).dat   = cdr.compileTimeParameteres;
dat(4).defaultNum = 0;

maxStrLen = 0;
for i = 1:length(dat)
  if ~isempty(dat(i).dat)
    maxStrLen = max([maxStrLen max(cellfun('length',{dat(i).dat(:).name}))]);
  end
end

for i = 1:length(dat)
  fprintf(f, '%s', prefix); fprintf(f, ' %s\n', dat(i).desc);
  dd = dat(i).dat;
  for j = 1:length(dd)
    fprintf(f, '%s', prefix); fprintf(f, '  %*s ', maxStrLen, dd(j).name);
    if i == 1 % Inputs
      sz = sprintf('%3i x %-3i', dd(j).size(1), dd(j).size(2));
    elseif i == 2 % Outputs
      sz = sprintf('%3s x %-3s', '?', '?');
    else
      eval(sprintf('varSet = ~isempty(cdr.%s);', dd(j).name));
      if j < dat(i).defaultNum || varSet
        eval(sprintf('[m,n] = size(cdr.%s);',dd(j).name));
        sz = sprintf('%3i x %-3i', m,n);
      else
        sz = sprintf('%9s','[SET ME]');
      end
    end
    fprintf(f, '%s', prefix); fprintf(f, '%s  %s\n', sz, dd(j).desc);
  end
  fprintf(f, '%s', prefix); fprintf(f, '\n');
end
