function gen_mv_exhaustive(A,~,~)
% Write out every non-zero multiplication
%
% Will be a stupid method for large matrices
%
%  y = A*x


pl('void mv_exhaustive(double y[%i], const double x[%i]) {', size(A,1), size(A,2))

A = full(A);

for r=1:size(A,1)
  p('y[%i] = 0.0', r-1)
  
  if size(A,2) > 50
    [~,ind,val] = find(A(r,:));
    if ~isempty(ind)
      p('+ %.20g*x[%i]', [val;ind-1]);
    end
  else
    
    for c=1:size(A,2)
      if abs(A(r,c)) > 1e-6
        if A(r,c) == 1
          p(' + x[%i]', c-1)
        elseif A(r,c) == -1
          p(' - x[%i]', c-1)
        elseif A(r,c) < 0
          p(' - %.20g*x[%i]', abs(A(r,c)), c-1)
        elseif A(r,c) > 0
          p(' + %.20g*x[%i]', A(r,c), c-1)
        end
      end
    end
  end
  pl(';')
end

pl('}')
