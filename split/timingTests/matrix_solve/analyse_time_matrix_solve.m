clear

T = readtable('data.dat');

A{1} = T(strcmp(T.Name,'Invert'),:);
A{2} = T(strcmp(T.Name,'LDL'),:);
A{3} = T(strcmp(T.Name,'Banded'),:);

symbols = {'o','x','+','s','*','d','v','^','<','>','p','h'};
colors  = get(gca,'ColorOrder');

figure(1);
clf
for i = 1:length(A)
  h(i) = semilogy(A{i}.Horizon, A{i}.Time,'marker',symbols{i},'color',colors(i,:));
  hold on
end
grid on

legend(h,{'Invert', 'LDL', 'Banded'},'location','NorthWest');
xlabel('Horizon')
ylabel('Time (sec)')
title('Linear system solve for n=10, m=3')
axis tight

prepFig
printFig('speed_horizon')


return


figure(1);
clf
for i = 1:length(A)
  dd = unique(A{i}.Density);
  for id = 1:length(dd)
    I = A{i}.Density == dd(id);
    loglog(A{i}.Rows(I), A{i}.Time(I),'marker',symbols{id},'color',colors(i,:));
    hold on
  end
end
grid on

for i = 1:length(A)
  h(i) = plot(10*[1;1],1e-6*[1;1],'color',colors(i,:),'linewidth',5);
end
dd = unique(A{1}.Density);
str = {};
for i = 1:length(dd)
  h(end+1) = plot(10*[1;1],1e-6*[1;1],'linestyle','none','marker',symbols{i},'color','k');
  str{i} = sprintf('Density %.2f', dd(i));
end
legend(h,{'Invert', 'LDL', str{:}},'location','NorthWest');
xlabel('Rows/Cols in A')
ylabel('Time (sec)')
title('Linear system solve product')
% str = sprintf('Density : ');
% title(str)
% title(sprintf('Computation time, Density %f', A{1}.Density(1)))
axis tight

prepFig
printFig('speed_density')

% return
% 
% figure(2);
% clf
% for i = 1:length(A)
%   loglog(A{i}.Rows, A{i}.Size/1e6);
%   hold on
% end
% grid on
% 
% legend({'Sparse', 'Dense', 'Exhaustive', 'BLAS'});
% xlabel('Rows/Cols in A')
% ylabel('Compiled size (MB)')
% title(sprintf('Program size, Density %f', A{1}.Density(1)))
% axis tight
% 
% prepFig
% printFig('size_density')
