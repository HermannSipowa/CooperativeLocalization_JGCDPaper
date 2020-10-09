close all
clear variables
clc
start_up
c1 = rgb('RosyBrown'); c2 = rgb('Black'); c3 = rgb('Lime');
c4 = rgb('Tomato'); c5 = rgb('DarkBlue'); c6 = rgb('DarkTurquoise');
ColorMatrix = [c1;c6;c3;c4;c5;c2];
Num_deputies = 5;

CentralizedFiterData = matfile('CentralizedData.mat');
Centralized_TracePos = CentralizedFiterData.PosData;

MonteCarloFilterData = matfile('MonteCarloTraceData.mat','Writable',true);
Agent1_ConstAnalysis = matfile('Agent1_ConstAnalysis.mat','Writable',true);
Agent2_ConstAnalysis = matfile('Agent2_ConstAnalysis.mat','Writable',true);
Agent3_ConstAnalysis = matfile('Agent3_ConstAnalysis.mat','Writable',true);
Agent4_ConstAnalysis = matfile('Agent4_ConstAnalysis.mat','Writable',true);
Agent5_ConstAnalysis = matfile('Agent5_ConstAnalysis.mat','Writable',true);
idxY = matfile('indexYCentralized.mat'); idxYCentralized = idxY.indexY;
idxY = matfile('indexYMonteCarlo.mat');  idxYMonteCarlo = idxY.indexY;
MCs = 724;
MC_intital = 1;
MC_end = MC_intital+MCs-1;
m   = size(MonteCarloFilterData.PosData_1,2);
n   = 6;
p   = 5;
jump = MCs*(0:5);
MCs_Trace(5).Data = zeros(MCs,m);
ConstAnalysis(5).Data = zeros(MCs,m);
for i = MC_intital:MC_end
    
    ii = i - MC_intital + 1;
    idx = ii + jump; % 1+n*(i-1):n*i;
    Posfield = strcat('PosData_',num2str(i));
    Constfield = strcat('Run_',num2str(i));
    TraceData = MonteCarloFilterData.(Posfield);
    for j = 1:p
        if j == 1
            MCs_Trace(j).Data(ii,:) = TraceData(1,:);
            ConstAnalysis(j).Data(idx,:) = Agent1_ConstAnalysis.(Constfield);
        elseif j == 2
            MCs_Trace(j).Data(ii,:) = TraceData(2,:);
            ConstAnalysis(j).Data(idx,:) = Agent2_ConstAnalysis.(Constfield);
        elseif j == 3
            MCs_Trace(j).Data(ii,:) = TraceData(3,:);
            ConstAnalysis(j).Data(idx,:) = Agent3_ConstAnalysis.(Constfield);
        elseif j == 4
            MCs_Trace(j).Data(ii,:) = TraceData(4,:);
            ConstAnalysis(j).Data(idx,:) = Agent4_ConstAnalysis.(Constfield);
        elseif j == 5
            MCs_Trace(j).Data(ii,:) = TraceData(5,:);
            ConstAnalysis(j).Data(idx,:) = Agent5_ConstAnalysis.(Constfield);
        end
    end
end


%% Plotting the trace of the position components
fh = figure; %('Renderer', 'painters', 'Position', [10 10 1000 700])
plt   = zeros(MCs, Num_deputies);
Label = {'Trace of P$_1$','Trace of P$_2$','Trace of P$_3$','Trace of P$_4$','Trace of P$_5$'};
for k = 1:Num_deputies
    subplot(3,2,k);
    plt(:,k) = plot(idxYMonteCarlo,MCs_Trace(k).Data,'Color',ColorMatrix(k,:));
    hold on
    plt2     = plot(idxYCentralized,Centralized_TracePos(k,:),'Color',c2);
    set(gca, 'YScale', 'log')
    ax = gca;
    ax.FontSize = 30;
    grid on
    xlabel('Period', 'FontSize', 30)
    ylabel(Label(k), 'FontSize', 30)
end
hL = legend([plt(end,1),plt(end,2),plt(end,3),plt(end,4),plt(end,5), plt2],...
    {'Agent1','Agent2','Agent3','Agent4','Agent5','Centralized Filter'});
hL.FontSize = 27;
newPosition = [0.68 0.17 0.1 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',2);


% set all units inside figure to normalized so that everything is scaling accordingly
 set(findall(fh,'Units','pixels'),'Units','normalized');
  % do not show figure on screen
%  set(fh, 'visible', 'off')
 % set figure units to pixels & adjust figure size
 fh.Units = 'pixels';
 fh.OuterPosition = [10 10 1000 700];
 % define resolution figure to be saved in dpi
 res = 500;
 % recalculate figure size to be saved
 set(fh,'PaperPositionMode','manual')
 fh.PaperUnits = 'inches';
 fh.PaperPosition = [0 0 7680 4320]/res;
 % save figure
 print(fh,'TraceMonteCarlo','-dpng',sprintf('-r%d',res))
 
%% Consistency Analysis
close all
Label = {'Agent1','Agent2','Agent3','Agent4','Agent5'};
fh2 = figure; %('Renderer', 'painters', 'Position', [10 10 1000 700]);
linesNum = 1:MCs;
for k = 1:Num_deputies
    plt = zeros(MCs);
    subplot(n/2,2,k)
    plt2 = yline(3,'--','LineWidth', 2, 'Color', 'r');
    yline(-3,'--','LineWidth', 2, 'Color', 'r');
    hold on
    for i = 1:n
        idx = linesNum + jump(i);
        plt(:,i) = plot(idxYMonteCarlo,ConstAnalysis(k).Data(idx,:),'Color',ColorMatrix(i,:));
    end
    ax = gca;
    ax.FontSize = 30;
    ylim([-10 10])
    xlabel('Period', 'FontSize', 30)
    ylabel(Label(k), 'FontSize', 30)
    hold off
    grid on
end
hL = legend([plt(end,1),plt(end,2),plt(end,3),plt2,plt(end,4),plt(end,5),plt(end,6)],...
    {'$\delta x$[km]', '$\delta y$[km]', '$\delta z$[km]','3$\sigma$ Covariance',...
    '$\delta \dot{x}$[m/s]', '$\delta \dot{y}$[m/s]', '$\delta \dot{z}$[m/s]'});
hL.FontSize = 27;
newPosition = [0.68 0.17 0.1 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',2);

% set all units inside figure to normalized so that everything is scaling accordingly
 set(findall(fh2,'Units','pixels'),'Units','normalized');
  % do not show figure on screen
 set(fh2, 'visible', 'off')
 % set figure units to pixels & adjust figure size
 fh2.Units = 'pixels';
 fh2.OuterPosition = [10 10 1000 700];
 % define resolution figure to be saved in dpi
 res = 500;
 % recalculate figure size to be saved
 set(fh2,'PaperPositionMode','manual')
 fh2.PaperUnits = 'inches';
 fh2.PaperPosition = [0 0 7680 4320]/res;
 % save figure
 print(fh2,'ConsistencyMonteCarlo','-dpng',sprintf('-r%d',res))

