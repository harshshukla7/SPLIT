function parseProb(cdr, prob)
dat = prob.dat;

if ~any(isfield(dat, {'pF','pB','pL'}))
  error('Can only generate solvers for paramertic problems')
end

% Fix matrix sizes
dat.f  = dat.f(:);
dat.pF = dat.pF';

% Add empty Argetric matrices if they're not present
if ~isfield(dat, 'pB')
  dat.pB = zeros(size(dat.A,1),0);
  dat.pF = zeros(size(dat.A,2),0);
end

% Build the prox functions
prox    = prob.prox;
nProx = length(prox);
L = []; l = []; pL = []; w = [];
proxInd = 0;
for i = 1:nProx
  L  = [L;prox(i).L];
  l  = [l;prox(i).l];
  proxInd(i+1) = size(L,1);
  if isfield(dat, 'pB')
    pL = [pL;full(prox(i).pL)];
  end
  dat.prox(i).type = prox(i).typeName;
  dat.prox(i).ind  = proxInd(i);
  dat.prox(i).len  = size(prox(i).L,1);
  dat.prox(i).dat  = prox(i).dat;
end
L = full(L); l = full(l);
if ~isfield(dat, 'pB')
  pL = zeros(size(L,1),0);
end

dat.L  = L;
dat.l  = l;
dat.pL = pL;
dat.proxInd = proxInd;

cdr.data = dat;
