clear all

d = splitData;

for i = 1:5
  if randn < 0
    d.add(sprintf('x%i', i), randn(ceil(rand*10),1));
  else
    d.add(sprintf('x%i', i), ceil(100*rand(ceil(rand*10),1)), 'int');
  end
end

d.writeFile;

