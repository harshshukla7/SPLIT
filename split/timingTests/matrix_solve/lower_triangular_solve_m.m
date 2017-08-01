function x = lower_triangular_solve_m(L, Lstr, s, p, b)



% % Store the strictly lower triangular matrix row-wise
% lvec = [];
% for i = 1:size(L,1)
%   lvec = [lvec L(i,1:i-1)];
% end

% Store the inverse of the diagonal
iDiagL = 1./diag(L);

Lt = tril(L,-1);
[nz_per_row,ind,vec] = splitData('test').writeSparseMatrix(Lt,'L')



i = 1;
for r = 1:size(L,1)
  x(r) = s(r) * b(p(r));
  for j = 1:nz_per_row(r)
    x(r) = x(r) - vec(i)*x(ind(i));
    i = i + 1;
  end
%   for c = 1:r-1
%     x(r) = x(r) - lvec(i)*x(c);
%     i = i + 1;
%   end
  x(r) = x(r) * iDiagL(r);
end


% i = 1;
% for r = 1:size(L,1)
%   x(r) = b(r);
%   for c = 1:r-1
%     x(r) = x(r) - lvec(i)*x(c);
%     i = i + 1;
%   end
%   x(r) = x(r) * iDiagL(r);
% end
