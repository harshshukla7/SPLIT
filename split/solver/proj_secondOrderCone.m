function x = proj_secondOrderCone(z)
%
% x = proj_secondOrderCone(z)
%
% Project the point z = [t;y] onto the second order cone
%
% x = argmin ||x - z||_2 s.t. ||y||_2 <= t
%

% Solution from Boyd's lecture notes
% http://www.seas.ucla.edu/~vandenbe/236C/lectures/proxop.pdf


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
