function sgn = sign(obj)
%
% Return the sign of the variable
%  
%  sgn = sign(x)
%
% sgn is a matrix of size x.size
% Returns NaN if the sign is unknown
%

sgn   = reshape(sign(obj.basis_(1,:)),obj.size_);
const = obj.isconstant;
nn    = NaN*ones(obj.size_); nn(const)=0;
sgn = sgn.*const + nn;
