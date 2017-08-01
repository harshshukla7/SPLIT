function solver = genMatlab(cdr)
% Generate a matlab function to solve the problem

fprintf('Generating matlab function...\n');

% Confirm that all required variables have been set
for i = cdr.numDefaultParams:length(cdr.inConstArgs)
  eval(sprintf('val = cdr.%s;', cdr.inConstArgs(i).name));
  if isempty(val)
    fprintf('Error: The parameter %s must be set before generating code\n\n', cdr.inConstArgs(i).name);
    return
  end
end
for i = 1:length(cdr.compileTimeParameteres)
  eval(sprintf('val = cdr.%s;', cdr.compileTimeParameteres(i).name));
  if isempty(val)
    fprintf('Error: The parameter %s must be set before generating code\n\n', cdr.compileTimeParameteres(i).name);
    return
  end
end

% Clear variables set during the generation (so they're not duplicated if
% genMatlab is called twice)
cdr.functions   = [];
cdr.allocVars   = [];
cdr.headerFuncs = [];

cdr.openGenFile(true);

template  = cdr.loadTemplate;

% Parse the file into tmp
ftmp = fopen('tmpx.m', 'w+');
assert(ftmp > 0, 'Could not open temporary file')
indIO     = find(~cellfun('isempty',regexpi(template,'\s*%+\s*>+\s*IO\s*>+.*')));
mode      = 0; % Copy-mode
for i = indIO+1:length(template)
  ln    = template{i};
  copy_line = ~isempty(regexpi(ln,'\s*%+\s*>+\s*COPY\s*>+.*'));
  run_line  = ~isempty(regexpi(ln,'\s*%+\s*>+\s*RUN\s*>+.*'));
  
  if     copy_line, mode = 0;
  elseif run_line,  mode = 1;
  else
    if mode == 0
      % Write the line into tmp as a print statement
      fprintf(ftmp, 'cdr.print(''%%s\\n'', ''%s'');\n', strrep(ln, '''',''''''));
    else
      % Write the line into tmp as a statement to exectue
      fprintf(ftmp, '%s\n', ln);
    end
  end
end
fclose(ftmp);

rehash % Tells matlab about the new tmp.m file
tmpx
delete('tmpx.m')

% Load the generated file into a text cell array
cdr.closeGenFile();

f = fopen(cdr.getFileName(),'r');
c = textscan(f, '%s', 'delimiter','\n','whitespace','');
fclose(f);

% fseek(cdr.fHdl, 0, 'bof')
% c = textscan(cdr.fHdl, '%s', 'delimiter','\n','whitespace','');
c = c{1};

% Find row where the >>> HEADER >>> line is
indHdr = find(~cellfun('isempty',regexpi(c,'\s*%+\s*>+\s*HEADER\s*>+.*')));
assert(length(indHdr) == 1, 'Invalid template file: Could not find >>>> HEADER >>>> row')

% Open yet again, clearing the file this time
cdr.openGenFile(false, 'w+');
fseek(cdr.fHdl, 0, 'bof');
cdr.genFunctionHeader;
for i = 1:indHdr-1,     cdr.print('%s\n', c{i}); end % Copy everything before the header
for i = 1:length(cdr.headerFuncs), feval(cdr.headerFuncs{i}); end % Process the header functions
for i = indHdr+1:length(c), cdr.print('%s\n', c{i}); end % Copy everything after the header

cdr.closeGenFile;




% Provide a function handle to call the function
str = sprintf('solver = @(%s)%s(%s,%s);',...
  cdr.csvList({cdr.inArgs.name}),...
  cdr.fname,...
  cdr.csvList({cdr.inArgs.name}),...
  cdr.csvList(cellfun(@(x)sprintf('cdr.%s',x), {cdr.inConstArgs.name},'UniformOutput',false)));
eval(str)
end
