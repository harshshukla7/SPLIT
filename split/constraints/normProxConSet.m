classdef normProxConSet < conSet
 % Prox-operator for a norm:
 % prox(v; lambda) = min_x norm(x, p) + (1/2*lambda)*||x-v||_2^2
 %                 = v - lambda*proj(v/lambda; p)
 
 properties
  type     = conSetTypes.normProx
  typeConj = conSetTypes.normBall
  p
  pDual
  cWeight
 end
 
 methods
  function con = normProxConSet(x, p)
   con.aff = x;
   con.p   = p;
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
   % dat.p - norm involved
   %
   
   prox = proxData@conSet(con);
   prox.dat.p = con.p;
   p = con.p;
   
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
   
   % specific value for the function proj_normBall ... 
   % it appears only in case of a norm in the objective
   
   prox.dat.c = prox.weight;
   
   prox.func = @(x, t) prox_norm(x, t, p);
   prox.conjFunc = @(x, c) proj_normBall(x, pDual, c);
   
   
   str = sprintf('prox_norm(%%s, %%s, %i)', con.p);
   prox.funcStr = @(x,rho) sprintf(str, x, rho);
   
   str = sprintf('proj_normBall(%%s, %i, %%s)', con.pDual);
   prox.conjStr = @(x,rho) sprintf(str, x, rho);
   
  end
 end
end