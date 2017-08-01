classdef normBallConSet < conSet
 % Constraint: norm-ball
 % norm(x, p) <= c
 
 properties
  type = conSetTypes.normBall
  typeConj = conSetTypes.normProx
  p
  pDual
  c
  sigma
 end
 
 methods
  function con = normBallConSet(x, p, c)
   con.aff = x;
   con.p   = p;
   con.c   = c;
   con.epi = x.getVars(true);
  end
  
  function prox = proxData(con)
   % prox = proxData(con)
   %
   % Return the data required to describe this prox operator
   %
   % Elements of struct:
   %  L, l - defining the affine portion of the constraint L*x+l in set
   %  dat  - set-specific data
   %  type - type of the constraint
   %
   % dat.c - level of the level-set
   %
   
   p = con.p; c = con.c; 
   prox = proxData@conSet(con);
   prox.dat.c = c;
   prox.dat.p = p;
   
   % I assign pDual here
   if p == 1
       pDual = Inf;
   elseif p == 2
       pDual = 2;
   elseif p == Inf || p == inf
       pDual = 1;
   else
       error('Unknown norm type');
   end
   prox.dat.pDual = pDual;
   con.pDual = pDual;

   prox.func = @(x,~) proj_normBall(x, p, c);
   prox.conjFunc = @(x, t) prox_norm(x, t, pDual);

   str = sprintf('proj_normBall(%%s, %i, %f)', con.p, con.c);
   prox.funcStr = @(x,~) sprintf(str, x);

   str = sprintf('proj_norm(%%s, %%s, %i)', con.pDual);
   prox.conjStr = @(x,rho) sprintf(str, x, rho);
  end
 end
end