classdef coderFunc_prox_conj < coderFunc
  
  methods
    function f = coderFunc_prox_conj(dat)
      % Create a function that evaluates the prox functions for the data
      %
      % y = prox(x)
      %
      
      f = f@coderFunc('void custom_prox(double y[nDual], double x[nDual])');
      
      for i = 1:length(dat.prox)
        prox = dat.prox(i);
        args = sprintf('y+%i, x+%i', prox.ind, prox.ind);
        %warning('Just created the coderFunc_prox_CP subclass from coderFunc_prox class. The only changed made is in swtich case (now case is swtiched by prox.typeConj.name). Please verify everything in coderFunc..._CP');
        switch prox.typeConj  %% This is the change
          case 'box',                 error ('Prox not implemented yet');
          case 'ellipse',             error ('Prox not implemented yet');
          case 'ellipseConj',         error ('Prox not implemented yet');
          case 'secondOrderCone'
            f.pl('proj_secondOrderCone(%s, %i);',args, prox.len)
          case 'secondOrderConeConj'
            f.pl('proj_secondOrderCone_conj(%s, %i);',args, prox.len)
          case 'nonPositive',
            f.pl('proj_negative(%s, %i);',args, prox.len)
          case 'nonNegative',
            f.pl('proj_positive(%s, %i);',args, prox.len)
          case 'normBall'
            switch prox.dat.p
              case 1,    funcName = 'proj_normball_one';
              case 2,    funcName = 'proj_normball_two';
              case inf,  funcName = 'proj_normball_inf';
              otherwise, error('Unkown norm-type')
            end
            f.pl('%s(%s, %g, %i);', funcName, args, prox.dat.c, prox.len)
          case 'normProx',
            switch prox.dat.p
              case 1,    funcName = 'prox_norm_one';
              case 2,    funcName = 'prox_norm_two';
              case inf,  funcName = 'prox_norm_inf';
              otherwise, error('Unkown norm-type')
            end
            f.pl('prox_var = %g*rho;',full(prox.dat.c))%
            
            f.pl('%s(%s, prox_var, %i);', funcName, args, prox.len)
            %f.pl('%s(%s, rho, %i);', funcName, args, prox.len)
          otherwise
            error('Unknown prox type')
        end
      end
    end
  end
end

%% To DO 

% if strcmp(char(prox(i).typeConj), 'normProx')
%     %          proxWeight(i) = rho*prec.P( proxInd{i}, proxInd{i}) * prox(i).dat.c;
%     proxWeight(i) = rho * prox(i).dat.c;
% end
% if strcmp(char(prox(i).typeConj),  'ellipseConj')
%     %          proxWeight(i) = rho*prec.P( proxInd{i}, proxInd{i}) * sqrt(prox(i).dat.c);
%     proxWeight(i) = rho * sqrt(prox(i).dat.c);
% end
