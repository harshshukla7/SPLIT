classdef nonnegConSet < conSet
 % Constraint: Positive orthant
 
 properties
  type = conSetTypes.nonNegative
  typeConj = conSetTypes.nonPositive
 end
 
 methods
  function con = nonnegConSet(x)
   con.aff = x;
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
   
   prox = proxData@conSet(con);
   prox.func = @(x,~) proj_positive(x);
   prox.conjFunc = @(x,~) proj_negative(x);
   
   prox.funcStr = @(x,~) sprintf('proj_positive(%s)', x);
   prox.conjStr = @(x,~) sprintf('proj_negative(%s)', x);   
  end
 end
end