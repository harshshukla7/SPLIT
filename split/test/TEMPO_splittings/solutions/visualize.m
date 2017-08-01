function [] = visualize(time, state_sim, l)

clf
whitebg([1.0 1.0 1.0])
set(gcf,'Color',[1 1 1])
hold on;

text(0.6,1.2,['current time: ' num2str(time(end)) 's'],'FontSize',15);
plot(state_sim(1,end),0,'ks','MarkerSize',30,'Linewidth',3);

theta = state_sim(2,end);
xB = state_sim(1,end)-l*sin(theta);
yB = l*cos(theta);

line([state_sim(1,end) xB], [0 yB], 'Linewidth',2);
plot(xB,yB,'ro','MarkerSize',25,'Linewidth',3);

grid on;
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);