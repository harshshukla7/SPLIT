classdef ellipseConSet < conSet
 % Constraint: ellipse
 %   x'*sQ'*sQ*x <= c
 
 properties
  type = conSetTypes.ellipse
  typeConj = conSetTypes.ellipseConj
  sQ
  c
 end
 
 methods
  function con = ellipseConSet(x, sQ, c)
   % ellipseConSet(x, sQ, c)
   %  x'*sQ'*sQ*x <= c
   con.aff = sQ*x;
   con.epi = x.getVars(true);
   
   con.sQ  = sQ;
   con.c   = c;
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
   
   prox = proxData@conSet(con);
   prox.dat.c = con.c;
   prox.dat.pDual = 2;
   
   pDual = prox.dat.pDual;
   sc = sqrt(con.c);
   
   prox.func = @(x,~) proj_quadball(x, sc);
   prox.conjFunc = @(x, t) prox_norm(x, t, pDual);
   
   prox.funcStr = @(x,~)   sprintf(sprintf('proj_quadball(%%s, %f)', full(sqrt(con.c))), x);
   prox.conjStr = @(x,rho) sprintf(sprintf('prox_norm(%%s, %%s, %i)', 2), x, rho);
  end
 end
end