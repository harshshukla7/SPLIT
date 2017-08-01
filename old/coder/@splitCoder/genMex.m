function genMex(cdr, cfg)
%----------------------------------------------------------------
% Call matlab coder to generate c-code for a mex file
%
%  genMex(cdr, cfg)
%
% cfg is a coder.config object (mex by default)
%----------------------------------------------------------------

% Generate the m-file first
cdr.genMatlab;

fprintf('Compiling mex-function...\n');

% Generate C-function from matlab file
% Rest of the arguments are constants
for i = 1:length(cdr.inArgs)
  args{i} = zeros(cdr.inArgs(i).size);
end
for i = 1:length(cdr.inConstArgs)
  eval(sprintf('args{end+1} = coder.Constant(cdr.%s);', cdr.inConstArgs(i).name));
end

% Add the include directory to the path
[pathstr,~,~]=fileparts(mfilename('fullpath'));
includePath = [pathstr filesep 'include'];

if nargin < 2
  cfg = coder.config('mex');
  cfg.ConstantInputs = 'Remove';
  cfg.IntegrityChecks = false;
  cfg.ResponsivenessChecks = false;
  cfg.DynamicMemoryAllocation = 'off';
end
str = ['codegen ' cdr.fname ' -args args -config cfg -I ' includePath ' splitTimer.h splitTimer.c;'];
eval(str)

