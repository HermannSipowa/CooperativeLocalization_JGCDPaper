close all
clc
c1 = rgb('RosyBrown'); c2 = rgb('Black'); c3 = rgb('Lime');
c4 = rgb('Tomato'); c5 = rgb('DarkBlue'); c6 = rgb('DarkTurquoise');
ColorMatrix = [c1;c6;c3;c4;c5;c2];
Num_deputies = 5;

CentralizedFiterData = matfile('CentralizedData.mat');
Centralized_TracePos = CentralizedFiterData.PosData;
Centralized_TraceVel = CentralizedFiterData.VelData;

MonteCarloFilterData = matfile('MonteCarloTraceData.mat');
Agent1_ConstAnalysis = matfile('Agent1_ConstAnalysis.mat');
Agent2_ConstAnalysis = matfile('Agent2_ConstAnalysis.mat');
Agent3_ConstAnalysis = matfile('Agent3_ConstAnalysis.mat');
Agent4_ConstAnalysis = matfile('Agent4_ConstAnalysis.mat');
Agent5_ConstAnalysis = matfile('Agent5_ConstAnalysis.mat');
idxY = matfile('indexYCentralized.mat'); idxYCentralized = idxY.indexY;
idxY = matfile('indexYMonteCarlo.mat');  idxYMonteCarlo = idxY.indexY;
MCs = 724;
m   = size(MonteCarloFilterData.PosData_1,2);
n   = 6;
p   = 5;
MCs_Trace(5).Data = zeros(MCs,m);
ConstAnalysis(5).Data = zeros(MCs,m);
for i = 1:MCs
    idx = 1+n*(i-1):n*i;
    pdx = 1+p*(i-1):p*i;
    Posfield = strcat('PosData_',num2str(i));
    Constfield = strcat('Run_',num2str(i));
    TraceData = MonteCarloFilterData.(Posfield);
    for j = 1:p
        if j == 1
            MCs_Trace(j).Data(i,:) = TraceData(1,:);
            ConstAnalysis(j).Data(idx,:) = Agent1_ConstAnalysis.(Constfield);
        elseif j == 2
            MCs_Trace(j).Data(i,:) = TraceData(2,:);
            ConstAnalysis(j).Data(idx,:) = Agent2_ConstAnalysis.(Constfield);
        elseif j == 3
            MCs_Trace(j).Data(i,:) = TraceData(3,:);
            ConstAnalysis(j).Data(idx,:) = Agent3_ConstAnalysis.(Constfield);
        elseif j == 4
            MCs_Trace(j).Data(i,:) = TraceData(4,:);
            ConstAnalysis(j).Data(idx,:) = Agent4_ConstAnalysis.(Constfield);
        elseif j == 5
            MCs_Trace(j).Data(i,:) = TraceData(5,:);
            ConstAnalysis(j).Data(idx,:) = Agent5_ConstAnalysis.(Constfield);
        end
    end
end


%% Plotting the trace of the position components
% 
% figure('Renderer', 'painters', 'Position', [10 10 1000 700])
% plt   = zeros(MCs, Num_deputies);
% for k = 1:Num_deputies
%     subplot(3,2,k);
%     plt(:,k) = plot(indexY,MCs_Trace(k).Data,'Color',ColorMatrix(k,:));
%     hold on
%     plt2     = plot(indexY,Centralized_TracePos(k,:),'Color',c2);
%     grid on
%     xlabel('Period')
%     ylabel(sprintf('Trace of P$_{%d}$', k))
%     set(gca, 'YScale', 'log')
% end
% hL = legend([plt(end,1),plt(end,2),plt(end,3),plt(end,4),plt(end,5), plt2],{'Agent1','Agent2','Agent3','Agent4','Agent5','Centralized Filter'});
% newPosition = [0.68 0.17 0.1 0.1];
% newUnits = 'normalized';
% set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',2);



%% Consistency Analysis

Label = {'Agent1','Agent2','Agent3','Agent4','Agent5'};
figure('Renderer', 'painters', 'Position', [10 10 1000 700])
for ii = 1:MCs
    idx = 1+n*(ii-1):n*ii;
    for k = 1:Num_deputies
        plt = zeros(6);
        subplot(n/2,2,k)
        plt2 = yline(3,'--','LineWidth', 2, 'Color', 'r');
        yline(-3,'--','LineWidth', 2, 'Color', 'r');
        hold on
        for i = 1:n
            plt(:,i) = plot(idxYMonteCarlo,ConstAnalysis(k).Data(idx(i),:),'Color',ColorMatrix(i,:));
        end
        xlabel('Period')
        ylabel(Label(k))
        hold off
        grid on
        clear StateError Sigma2
    end
end
hL = legend([plt(end,1),plt(end,2),plt(end,3),plt2,plt(end,4),plt(end,5),plt(end,6)],...
    {'$\delta x$[km]', '$\delta y$[km]', '$\delta z$[km]','3$\sigma$ Covariance',...
    '$\delta \dot{x}$[m/s]', '$\delta \dot{y}$[m/s]', '$\delta \dot{z}$[m/s]'});
newPosition = [0.68 0.17 0.1 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',2);


