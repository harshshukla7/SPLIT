function func_mldivide(dat, func_name, A, varargin)
%
% func_mldivide(dat, func_name, A, param1, val1, etc)
%
% Generate code to solve a system of equations
%
% y = A\x
%
% Function prototype:
%   void func_name(y[n], x[n])
%
% Arguments:
%  func_name  : string
%  A          : Matrix data
%
% Optional:
%  Astr       : Name to use when storing the matrix [otherwise random]
%  method     : One of 'ldl', 'invert'
%               If not specified, method will be chosen based on the
%               structure of A
%  desc       : Text description of the function
%  zero_tol   : Value below which a number is considered zero [1e-10]
%               Used to sparsify the invert method
%

p = inputParser;
addRequired(p, 'func_name',  @ischar);
addRequired(p, 'A',          @isnumeric);

addParameter(p, 'Astr',       dat.uniqueVarName, @ischar);
addParameter(p, 'method', 'auto', ...
  @(x) any(validatestring(x,{'auto','invert','ldl'})));
addParameter(p, 'desc', '', @ischar);
addParameter(p, 'zero_tol', 1e-10, @isnumeric);

% Parse and copy the results to the workspace
parse(p, func_name, A, varargin{:});
cellfun(@(f) evalin('caller',[f ' = p.Results.' f ';']), fieldnames(p.Results))

if strcmpi(method, 'auto')
  % Try and guess the right method
  warning('Harsh / Giorgos : Put code here to guess which approach should be used by analysing the matrix A and the timing plots')
  method = 'invert';
end

% Write the function prototype
[m,n] = size(A);
assert(m==n, 'func_mldivide currently only solves square matrices')
f = coderFunc('void %s(REAL y[%i], const REAL x[%i])', func_name, n, n);

% Write the code to actually solve the problem
switch method
  %%
  case 'invert'
    % Invert the matrix and solve the problem via matrix-vector
    % multiplication
    invA = inv(full(A));
    invA(abs(invA) < zero_tol) = 0;
    fInv = dat.func_times([func_name '_inv'], invA);
    
    % Copy to the body of the function into ours
    f.p(fInv.body);
    
    %%
  case 'ldl'
    % Solve the symmetric system of equations via LDL pre-factorization
    assert(norm(full(A - A')) == 0, 'Matrix must be symmetric to use LDL factorization')
    
    % Pre-factor the matrix A
    A = sparse(A);
    [L,D,p,S] = ldl(A,0.01,'vector');
    s = diag(S);
    s = s(p);
    
    % Create all the functions that we will need:
    
    % y  = L \ (s .* b(p));
    dat.functions(end+1) = ...
      lower_triangular_solve(dat, L, [Astr '_L'], s, p, [func_name '_solve1']);
    
    % z  = inv(D) * y;
    dat.functions(end+1) = ...
      dat.func_times([func_name '_solve2'], inv(D), 'method', 'exhaustive');
    
    % q  = L' \ z;
    % x(p,1)  = s .* q;
    dat.functions(end+1) = ...
      upper_triangular_solve(dat, L', [Astr '_LT'], s, p, [func_name '_solve3']);
    
    % Create the function that ties everything together
    f = coderFunc('void %s(double _x[%i], double _b[%i])', func_name, length(s), length(s));
    f.pl('  %s_solve1(_x,  _b);', func_name)
    f.pl('  %s_solve2(_b,  _x);', func_name)
    f.pl('  %s_solve3(_x,  _b);', func_name)
end

dat.functions(end+1) = f;
end


function f = lower_triangular_solve(dat, L, Lstr, s, p, funcname)
% funcname = forward_solve(dat, L, Lstr, s, p)
%
% Generate sparse code for forward solve of
%   L*x = s .* b(p)
%
% L    : lower triangular matrix
% Lstr : name of the matrix L (for dense generation)
% p    : permutation vector
% s    : scaling vector
%
% funcname(double *x, double *b)
%
% NOTE : b is changed on exit!
%

% Store the inverse of the diagonal
iLDiag = sprintf('%s_iDiag', Lstr);
dat.add_var(iLDiag, 1./diag(L));

% Store the lower-triangular portion of L in sparse form
Ltri = sprintf('%s_lower_tri', Lstr);
dat.add_var(Ltri, sparse(tril(L,-1)));

% % Store the permutation vector and scaling
pName = sprintf('%s_p', Lstr);
sName = sprintf('%s_s', Lstr);
dat.add_var(pName, p-1, 'type', 'int', 'storage_method', 'dense');
dat.add_var(sName, s,   'type', 'real', 'storage_method', 'dense');

f = coderFunc('void %s(double _x[%i], double _b[%i])', funcname, size(L,2), size(L,1));
f.pl('int i = 0;')
f.pl('for(int r=0; r<%i; r++) {', size(L,1));
f.pl('  _x[r] = %s[r] * _b[%s[r]];', sName, pName);
f.pl('  for(int j=0; j<%s_nz_per_row[r]; j++) {', Ltri);
f.pl('    _x[r] -= %s_dat[i] * _x[%s_ind[i]];', Ltri, Ltri);
f.pl('    i++;');
f.pl('  }');
f.pl('  _x[r] *= %s[r];', iLDiag);
f.pl('}');
end

function f = upper_triangular_solve(dat, U, Ustr, s, p, funcname)
% funcname = upper_triangular_solve(dat, U, Ustr, s, p, funcname)
%
% Generate sparse code for the solve of
%   x(p) = s .* (inv(U) * z)
%
% U    : upper triangular matrix
% Ustr : name of the matrix U
% p    : permutation vector
% s    : scaling vector
%
% funcname(double *x, double *z)
%
% NOTE : z is modified after call

%       q  = L' \ z;
%       x(p,1)  = s .* q;

% Store the inverse of the diagonal
iUDiag = sprintf('%s_iDiag', Ustr);
dat.add_var(iUDiag, 1./diag(U));

% Store the upper-triangular portion of U in sparse form
Utri = sprintf('%s_upper_tri', Ustr);
dat.add_var(Utri, sparse(triu(U,1)));

% Store the scaling vector
sName = sprintf('%s_s', Ustr);
dat.add_var(sName, full(s));

f = coderFunc('void %s(double _x[%i], double _z[%i])', funcname, size(U,2), size(U,1))
f.pl('int i = %s_ind_len-1;', Utri)
f.pl('for(int r=%i; r>=0; r--) {', size(U,1)-1);
f.pl('  _x[r] = _z[r];');
f.pl('  for(int j=0; j<%s_nz_per_row[r]; j++) {', Utri);
f.pl('    _x[r] -= %s_dat[i] * _x[%s_ind[i]];', Utri, Utri);
f.pl('    i--;');
f.pl('  }');
f.pl('}');

% Scale and permute
for i = 1:size(U,1)
  f.pl('  _z[%i] = %s[%i] * _x[%i];', p(i)-1, sName, i-1, i-1);
end

% Copy to output
f.pl('  memcpy(_x, _z, sizeof(double)*%i);', size(U,1));
end
