classdef conSet
% 
% conSet
%
% Represents a prox function
%
%  (Poorly named - origionally represented set constraints, but now also
%  includes all types of prox operators)
%
 
 properties
  epi  = []                   % epigraph variables that this constraint needs
  aff  = []                   % Affine term aff in Set
  weight = 1                  % Weighting multiplier for this function (doesn't matter if it is an indicator function)
 end

 properties (Abstract)
  type                        % Type of the proximal operator
  typeConj
 end 
 
 methods
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
   
   prox.type = con.type;
   prox.typeName = con.type.name;
   prox.typeConj = con.typeConj;
   
   B = con.aff.basis;
   prox.L   = B(2:end,:)';
   prox.l   = B(1,:)';
   prox.dat = [];
   prox.func = @(~,~) error('Prox function not defined for this constraint');
   prox.conjFunc = @(~,~) error('Conjugate not defined for this constraint');
   prox.weight = con.weight;
   
   % Strings to evaluate the prox and conjugate functions
   prox.funcStr = @(x,rho) error('Prox function not defined for this contraint');
   prox.conjStr = @(x,rho) error('Conjugate function not defined for this contraint');
  end
 end
end