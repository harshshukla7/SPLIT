classdef boxConSet < conSet
 % Constraint: norm-ball
 % norm(x, p) <= c
 
 properties
  type = conSetTypes.box
  lb % lower bound
  ub % upper bound
 end
 
 methods
  function con = boxConSet(x, lb, ub)
   con.aff  = x;
   con.lb   = lb;
   con.ub   = ub;
   con.epi  = x.getVars(true);
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
   
   lb = con.lb; ub = con.ub;
   prox = proxData@conSet(con);
   prox.dat.lb = lb;
   prox.dat.ub = ub;
   
   prox.func = @(x,~) proj_box(x, lb, ub);
   
   % Return a string that evaluates the proximal function
   lbStr = ['[' sprintf('%f;', con.lb) ']'];
   ubStr = ['[' sprintf('%f;', con.ub) ']'];
   str = sprintf('proj_box(%%s, %s, %s)', lbStr, ubStr);
   prox.funcStr = @(x,~) sprintf(str, x);
  end
 end
end