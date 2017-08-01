function tf = isaffine(x)
%
% tf = isaffine(x)
%
% Returns a logical matrix of size x
%  tf(i,j) is true is the function x(i,j) is affine and false otherwise
%

B = x.basis ~= 0;
B = B(2:end,:);

% Get all function types
[~,types] = splitProb.getFunctions;

tf = zeros(size(B,2),1) == 0;
for i = 1:size(B,2)
 if any(types(B(:,i)) ~= splitFuncTypes.affine)
  tf(i) = false;
 end
end

tf = reshape(tf,size(x,1),size(x,2));
