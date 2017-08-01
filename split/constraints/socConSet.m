classdef socConSet < conSet
 % Constraint: second-order cone
 % norm(x, 2) <= a'*x + c
 
 properties
  type = conSetTypes.lorentz
  typeConj = conSetTypes.lorentzConj
  x
  t
 end
 
 methods
  function con = socConSet(x, t)
   % socConSet(x, t)
   %  norm(x, 2) <= t
   con.aff = [t;x];
   con.x   = x;
   con.t   = t;
   con.epi = unique([x.getVars(true);t.getVars(true)]);
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
   
   prox = proxData@conSet(con);
   prox.func = @(x,~) proj_secondOrderCone(x);
   prox.conjFunc = @(x,~) proj_secondOrderConeConj(x);
   
   prox.funcStr = @(x,~) sprintf('proj_secondOrderCone(%s)', x);
   prox.conjStr = @(x,~) sprintf('proj_secondOrderConeConj(%s)', x);   
  end
 end
end