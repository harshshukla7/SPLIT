classdef coderFunc_prox < coderFunc
  
  methods
    function f = coderFunc_prox(dat)
      % Create a function that evaluates the prox functions for the data
      %
      % y = prox(x)
      %
      
      f = f@coderFunc('void custom_prox(real y[nDual], real x[nDual])');
      
      for i = 1:length(dat.prox)
        prox = dat.prox(i);
        args = sprintf('y+%i, x+%i', prox.ind, prox.ind);
        switch prox.type
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
            %f.pl('prox_var = %g/rho;',prox.dat.c/rho)
            f.pl('%s(%s, %g, %i);', funcName, args, prox.dat.c, prox.len)
          case 'normProx',
            switch prox.dat.p
              case 1,    funcName = 'prox_norm_one';
              case 2,    funcName = 'prox_norm_two';
              case inf,  funcName = 'prox_norm_inf';
              otherwise, error('Unkown norm-type')
            end
            f.pl('prox_var = %g/rho;',full(prox.dat.c))
            f.pl('%s(%s, prox_var, %i);', funcName, args, prox.len)
            %f.pl('%s(%s, rho, %i);', funcName, args, prox.len)
          otherwise
            error('Unknown prox type')
        end
      end
    end
  end
end
