function plotTACfig( res, parameter)
%% PLOTTACFIG Plots the figures for the TAC paper

M = res.sdpSol.M;
d = res.sdpSol.d;
Lambda = res.sdpSol.Lambda;

N = parameter.N;
model = res.model;
nu = size(model.ssM.Bu,2);
nd = size(model.ssM.Bd,2);
nr = size(model.outputs.D,1);
nrTot = nr*sum(parameter.active);
if strcmpi(parameter.omegaType,'slack')
    ns = nr;
else
    ns = 0;
end
nsTot = ns*sum(parameter.active);
% overall slightly confusing in terms of the different things computed 
% (input, slack, modified uncertainty with or without weather...)
ny = size(model.ssM.C,1);
if ~parameter.weather.uncertain
    F = res.Xi.F{1};
    f = res.Xi.f{1};
    unom = d(nu*N+nsTot+1:end)';
else
    F = res.Xi.F{end}(:,nd*N+1:end);
    f = res.Xi.f{end};
    unom = d((nu+nd)*N+nsTot+1:end)';
end
%% Generate trajectories for disturbance
Tot = 5;
switch parameter.ref.type
    case 'battery'
        nTraj = 0;
        trajUnc = zeros(nrTot,Tot);
        while nTraj<=Tot
            tr = 2*(rand(nrTot,1)-0.5);
            if (F*tr + f >=0)
%                 nTraj;
                nTraj = nTraj + 1;
                trajUnc(:,nTraj) = tr;
            end
        end
    case 'box'
        trajUnc = 2*(rand(nrTot,Tot)-0.5);
    case 'ellipsoid'
        % unit ball
        trajUnc = sampleUnitBall(Tot,nrTot);
end

if parameter.weather.uncertain
    % Generate trajectories for weather
    nominal = res.model.dist(:,1:N);
    
    %box
%     trajWeath = (1 + 0.1*(rand(nd,N,Tot)-0.5)).*nominal(:,:,ones(Tot,1));
    
    % In the direct product of the three ellipsoids
    normWeath = sampleUnitBall(3*Tot,N)';
    normWeath = reshape(normWeath,3,N,Tot);
    maxNom = max(nominal,[],2);
    trajWeath = nominal(:,:,ones(Tot,1)) + 0.1*maxNom(:,ones(N,1),ones(Tot,1)).*normWeath;
else
    trajWeath = res.model.dist(:,1:N,ones(Tot,1));
end

%% color definition
fancyred=[230	69	69	]/255;
fancyblue=[37 178 207]/255;
fancygreen=[77 200 161]/255;
fancyyellow=[232 205 35]/255;
orange=[255 175 42]/255;
darkorange=[200 100 50]/255;
fancycyan=[0.2 0.7 0.9];
orchid=[220 100 220]/255;

figure(1)
clf

%Activation times
act = find(parameter.active);

for i = 1:Tot
    weather = trajWeath(:,:,i);
    if parameter.weather.uncertain
        input = M*[weather(:); trajUnc(:,i)] + d;
    else
        input = M*[trajUnc(:,i)] + d;
    end
    u = reshape(input(1:nu*N),nu,[]);       %Inputs to the building
%     xi = reshape(input(nu*N+1:end),nr+nd,[]);   %Disturbances (reference)
    s = reshape(input(nu*N+1:nu*N+nsTot),ns,[]);
    w = input(nu*N+nsTot+1:end);    % reshaped disturbance
    if parameter.weather.uncertain
        dw = reshape(w(1:nd*N),nd,[]);
        r = reshape(w(nd*N+1:end),nr,[]);
    else
        r = reshape(w,nr,[]);
    end
    x = zeros(length(model.x0),N+1);        %States of the building
    x(:,1) = model.x0;
    y = zeros(ny,N);                        %Outputs of the building
    p = zeros(nr,N);                        %Total consumption
    for t =1:N
        x(:,t+1) = model.ssM.A*x(:,t)+model.ssM.Bu*u(:,t) + model.ssM.Bd*weather(:,t);
        y(:,t) = model.ssM.C*x(:,t);
        p(:,t) = model.outputs.D*u(:,t);
    end
    uACS = r-unom;                           %battery requirement
    s = cumsum(uACS);                       % state of the battery
    e = p(parameter.active)-r;                                   % error of the tracking
    %% Plotting
    %
    timeScale = 1:N;
    
    subplot(6,1,1);
    hold on;
    plotY = plot(y','color',[0.5 0.5 0.5],'linewidth',1.5);
    plotMeanY = plot(mean(y),'color',fancyred,'linewidth',1.5);
    %
    subplot(6,1,2);
    hold on;
    plotR = stairs(p,'color',orchid,'LineWidth',1.5);
    %
    subplot(6,1,3)
    hold on
    plotACS = stairs([act act(end)+1],[uACS uACS(end)],'color',fancyblue,'LineWidth',1.5);
    %
    subplot(6,1,4)
    hold on;
    plotS = stairs([act act(end)+1],[s s(end)],'color',fancygreen,'LineWidth',1.5);
    
    subplot(6,1,5)
    hold on;
    plotE = stairs([act act(end)+1],[e e(end)],'color',darkorange,'LineWidth',1.5);
    
    if i == 1
        subplot(6,1,6)
        [hAx,hLine1,hLine2] = plotyy(timeScale,weather(1,:),timeScale,weather(2,:));
        hold(hAx(1),'on')
        hold(hAx(2),'on')
        hLine1.Color = fancycyan;
        hLine2.Color = orange;
    else
%         hold(hAx(1),'on')
%         hold(hAx(2),'on')
        plotTmp = plot(hAx(1),timeScale,weather(1,:),'color',fancycyan);
        plotSun = plot(hAx(2),timeScale,weather(2,:),'color',orange);
    end
end

%% Prepare plots: add legends, limits...

subplot(6,1,1);
hold on
plotMaxY = plot([1 N+1],[ones(2,1)*max(model.constraints.ymax) ones(2,1)*min(model.constraints.ymin)],'k--');
h = rectangle('Position',[act(1), -100, act(end)-act(1)+1, 200],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
uistack(h,'bottom');
xlim([1 N+1]);
ylim([min((model.constraints.ymin))-0.5 max((model.constraints.ymax))+3]);
ylabel('Temperature [C]')
% xlabel('Time [h]')
legend([plotMeanY,plotY(1),plotMaxY(1)],{'Mean Tmp','Zone Tmp','Comfort'},...
        'Orientation','HORIZONTAL','Location','NorthEast');
legend('boxoff')

subplot(6,1,2)
hold on
h = rectangle('Position',[act(1), -100, act(end)-act(1)+1, 200],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
uistack(h,'bottom');
ylim([-42 10])
xlim([1 N+1])
legend([plotR],'Agg. input',...
    'Orientation','HORIZONTAL','Location','North');
legend('boxoff')
ylabel('Power (kW)')
% xlabel('Time [h]')
    
subplot(6,1,3)
hold on
h = rectangle('Position',[act(1), -100, act(end)-act(1)+1, 200],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
uistack(h,'bottom');
plotMaxPow = plot([1 N+1],[ones(2,1)*Lambda(1) -ones(2,1)*Lambda(1)],'k--');
ylim([-Lambda(1)-1 +Lambda(1)+4.5])
xlim([1 N+1])
legend([plotACS,plotMaxPow(1)],'Reference signal','Reference limits',...
    'Orientation','HORIZONTAL','Location','North');
legend('boxoff')
ylabel('Power (kW)')
% xlabel('Time [h]')


subplot(6,1,4)
hold on;
plotMaxS = plot([1 N+1],[2.8*ones(2,1)*Lambda(1) -2.8*ones(2,1)*Lambda(1)],'k--');
h = rectangle('Position',[act(1), -100, act(end)-act(1)+1, 200],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
uistack(h,'bottom');
ylim([-Lambda(1)*2.8-1 Lambda(1)*2.8+10])
xlim([1 N+1])
legend([plotS,plotMaxS(1)],'Reference integral','Capacity limits',...
    'Orientation','HORIZONTAL','Location','North');
legend('boxoff')
ylabel('Stored Energy (kWh)')
% xlabel('Time [h]')

subplot(6,1,5)
hold on;
plotMaxE = plot([1 N+1],[0.1*ones(2,1)*Lambda(1) -0.1*ones(2,1)*Lambda(1)],'k--');
h = rectangle('Position',[act(1), -100, act(end)-act(1)+1, 200],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
uistack(h,'bottom');
ylim([-Lambda(1)*0.1-0.2 Lambda(1)*0.1+0.6])
xlim([1 N+1])
legend([plotE,plotMaxE(1)],'Tracking Error','Max. error',...
    'Orientation','HORIZONTAL','Location','North');
legend('boxoff')
ylabel('Tracking error(kW)')
% xlabel('Time [h]')

set(hAx,{'ycolor'},{fancycyan;orange})
xlim(hAx(1),[1 N+1])
xlim(hAx(2),[1 N+1])
ylim(hAx(1),[10 30])
ylim(hAx(2),[-20 450])
set(hAx(2),'YTick',[0 200 400])
ylabel(hAx(1),'Outdoor Tmp. (C)')
ylabel(hAx(2),'Sun radiation(W/m^2)')
xlabel(hAx(1),'Time[h]')

cleanfigure; 
matlab2tikz('timePlotTAC.tex', 'height', '\figureheight', 'width', '\figurewidth')
    

end

