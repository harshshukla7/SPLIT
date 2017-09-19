classdef parameter < splitvar
  % Wrapper of the splitvar class that registers the splitvar as a
  % parameter with the splitProb
  %
  % Use the set function to set the current value of the parameter and then
  % the splitProb.genParamProb to evaluate
  
  methods
    function p = parameter(n, m)
      if nargin < 2, m = 1; end
      if nargin < 1, n = 1; end
      
      % Initialize splivar superclass
      p@splitvar(n,m);
      
      % Register the variable as a parameter
      splitProb.setVarToParameter(p.varNum_(:));
    end
  end
end
