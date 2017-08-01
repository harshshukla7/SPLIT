function testProx()

N = 100;

i = 1;
dat(i).func = @(z,c) prox_norm(z,c,1);
dat(i).name = 'prox_norm_one';
i = i + 1;
dat(i).func = @(z,c) prox_norm(z,c,2);
dat(i).name = 'prox_norm_two';
i = i + 1;
dat(i).func = @(z,c) proj_normBall(z,2,c);
dat(i).name = 'prox_normball_two';
i = i + 1;
dat(i).func = @(z,c) proj_normBall(z,inf,c);
dat(i).name = 'prox_normball_two';
i = i + 1;
dat(i).func = @(z,c) proj_secondOrderCone(z);
dat(i).name = 'proj_secondOrderCone';
i = i + 1;
dat(i).func = @(z,c) proj_negative(z);
dat(i).name = 'proj_negative';
i = i + 1;
dat(i).func = @(z,c) proj_positive(z);
dat(i).name = 'proj_positive';

for func = 1:length(dat)
  fprintf('Testing %-20s...', dat(func).name);
  for i = 1:N
    [y,c] = getRandVector();
    
    mm = dat(func).func(y,c);
    mc = test_matrix(y,c,func);
    
    err(i) = norm(mm - mc);
  end
  fprintf('  Max error = %f\n', max(err));
end
end

function [y,c] = getRandVector()
  MAXSIZE = 1000;
  MAXC = 10;

  m = ceil(rand*MAXSIZE);
  y = randn(m,1);
  c = rand*MAXC;
end
