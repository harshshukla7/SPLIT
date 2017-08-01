classdef conSetTypes
 properties
  name % Name that will be exported to describe this constraint
  num  % Unique number that will be used to describe this constraint
 end
 
 enumeration
  box          (1, 'box')
  ellipse      (2, 'ellipse')
  ellipseConj  (3, 'ellipseConj')
  lorentz      (4, 'secondOrderCone')
  lorentzConj  (5, 'secondOrderConeConj')
  nonPositive  (6, 'nonPositive')
  nonNegative  (7, 'nonNegative')
  normBall     (8, 'normBall')
  normProx     (9, 'normProx')
  noConstraint (10, 'ERROR - THIS SHOULD NEVER APPEAR')
 end

 methods
  function obj = conSetTypes(num, name)
   obj.num  = num;
   obj.name = name;
  end
 end
 
end