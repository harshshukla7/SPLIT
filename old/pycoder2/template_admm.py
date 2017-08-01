import splitCoder as split
import numpy as np
from   numpy.linalg import inv

nDual,nPrimal = split.dat['L'].shape
nEqCon        = split.dat['A'].shape[0]
nParam        = split.args['par']

# Arguments to the function
for key,val in split.args.iteritems():
  locals()[key] = split.Arg(key, val)

# Variable definitions
x        = split.Var('x',        nPrimal+nEqCon) # Extra rows are working space for solving KKT system
y        = split.Var('y',        nDual)
lam      = split.Var('lam',      nDual)
prev_lam = split.Var('prev_lam', nDual)
prev_y   = split.Var('prev_y',   nDual)

rDual    = split.Var('rDual',   dtype=float);
rPrimal  = split.Var('rPrimal', dtype=float);
itr      = split.Var('rDual',   dtype=int);
i        = split.Var('rPrimal', dtype=int);

# Load variables to compute parametric values
par_vars = (('l','pL'),('f','pF'),('b','pB'))
for var,parVar in par_vars:
  if split.dat[parVar].size == 0: # l = l_const
    locals()[var] = split.Const(var, split.dat[var])
  else:  # l = pL*par + l_const
    locals()[var]              = split.Var(var, nDual) # Variable l
    locals()[parVar]           = split.Const(parVar)   # Constant pL
    locals()['%s_const' % var] = split.Const(var)      # Constant l_const

# Pre-compute the inverse of the KKT matrix
# NOTE : This assumes that rho is constant
rho = split.constants['rho']
Q   = split.dat['Q']
L   = split.dat['L']
A   = split.dat['A']

# K = [Q+rho*L'*L A'; A zeros(size(A,1))];
K = split.concat([[Q+rho*np.dot(L.T,L), A.T], [A, np.zeros([A.shape[0]]*2)]])
invK = inv(K)
invK[np.abs(invK) < 1e-6] = 0

invK     = split.Const('invK', invK)

# Function header

split.gen(l, 'pL*par + l_const')

// Set kktRHS[nDual+1:end] = pB*par + b
copy_vector(kktRHS+nDual, b, nEqCon);

  for(itr=0; itr<MAXITR; itr++):
    prev_lam == lam
    prev_y   == y

    # Solve the KKT system
    # Solve using dense mat-vec product, since we pre-inverted the KKT matrix
    tmp2 == sum([-l, y, -lam])
    tmp  == L.T * tmp2
    tmp  == rho * tmp - f
    x    == invK * tmp

    split.gen('x = invK * (rho * (L.T * sum([-l, y, -lam]) - f)')

    workDual = lam + L*x + l

    # Evaluate prox functions 
    y = prox(workDual, rho)

    # Dual update
    lam = workDual - y;

    # Check convergence
    if (itr % ITR_PER_CONV_TEST == 0):
      rDual = sum_squared(lam - prev_lambda)

      if (rDual < DUALTOL) 
        rPrimal = sum_squared(rho*L*(y - prev_y))

        if (rPrimal < PRIMALTOL):
          break

  // Copy to solution structure
  copy_vector(sol->primal, x, nPrimal);
  copy_vector(sol->dual, lambda, nDual);
  sol->itr = itr;
  sol->rDual = rDual;
  sol->rPrimal = rPrimal;
}






# Initialize all variables to zero
# Define function prototype
