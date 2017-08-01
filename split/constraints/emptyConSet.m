classdef emptyConSet < conSet
 % Constraint: empty / no constraint
 % Used to pass epigraph constraints back to splitProb
 
 properties
  type = conSetTypes.noConstraint
 end
 
 methods
  function con = emptyConSet(x)
   con.epi   = x.getVars(true);
  end
 end
end