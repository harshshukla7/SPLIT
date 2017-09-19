function x = proj_socpconj1(z)
%
% x = proj_secondOrderConeConj(z)
%
% Project the point z = [t;y] onto the second order cone
%
% x = argmin ||x - z||_2 s.t. ||y||_2 <= t
%

% modify the input vector
z = [-z(1); z(2:end)];

t = z(1);
y = z(2:end);
ny = norm(y);

if ny <= t
  x = z;
elseif ny <= -t
  x = zeros(length(z), 1);
else
  x = [ny;y] * (t+ny)/(2*ny);
end

% modify the output vector
x = [x(2:end,1); -x(1,1)];