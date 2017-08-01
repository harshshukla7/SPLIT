function A = createbandOPT(n,b)
% b is the half of the bandwidth i.e. superdiagnoal or subdiagonal nozeros
row = [1:n];
col = [1:n];
for ii = 1:b
    le = length(row)+1;
    row(1,le:le+n-ii-1) = ii+1:n;
    col(1,le:le+n-ii-1) = 1:n-ii;
end
for ii = 1:b
    le = length(row)+1;
    row(1,le:le+n-ii-1) = 1:n-ii;
    col(1,le:le+n-ii-1) = ii+1:n;
end
A = sparse(row,col,1);