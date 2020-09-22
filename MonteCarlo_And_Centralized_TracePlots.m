close all
clc
c1 = rgb('RosyBrown'); c2 = rgb('Black'); c3 = rgb('Lime');
c4 = rgb('Tomato'); c5 = rgb('DarkBlue'); c6 = rgb('DarkTurquoise');
ColorMatrix = [c1;c6;c3;c4;c5;c2];
Label = {'Agent1','Agent2','Agent3','Agent4','Agent5'};
Num_deputies = 5;

%% Plotting the trace of the position components
MonteCarlo_TracePos(1).Data  = dlmread('Agent1_TracePosCovariance_MonteCarlo.dat');
MonteCarlo_TracePos(2).Data  = dlmread('Agent2_TracePosCovariance_MonteCarlo.dat');
MonteCarlo_TracePos(3).Data  = dlmread('Agent3_TracePosCovariance_MonteCarlo.dat');
MonteCarlo_TracePos(4).Data  = dlmread('Agent4_TracePosCovariance_MonteCarlo.dat');
MonteCarlo_TracePos(5).Data  = dlmread('Agent5_TracePosCovariance_MonteCarlo.dat');

Centralized_TracePos(1).Data  = dlmread('Agent1_TracePosCovariance_Centralized.dat');
Centralized_TracePos(2).Data  = dlmread('Agent2_TracePosCovariance_Centralized.dat');
Centralized_TracePos(3).Data  = dlmread('Agent3_TracePosCovariance_Centralized.dat');
Centralized_TracePos(4).Data  = dlmread('Agent4_TracePosCovariance_Centralized.dat');
Centralized_TracePos(5).Data  = dlmread('Agent5_TracePosCovariance_Centralized.dat');
figure('Renderer', 'painters', 'Position', [10 10 1000 700])

plt   = zeros(size(MonteCarlo_TracePos(1).Data,1), Num_deputies);
for k = 1:Num_deputies
    subplot(3,2,k);
    plt(:,k) = plot(indexY,MonteCarlo_TracePos(k).Data,'Color',ColorMatrix(k,:));
    hold on
    plt2     = plot(indexY,Centralized_TracePos(k).Data,'Color',c2);
    grid on
    xlabel('Period')
    ylabel(sprintf('Trace of P$_{%d}$', k))
    set(gca, 'YScale', 'log')
end
hL = legend([plt(end,1),plt(end,2),plt(end,3),plt(end,4),plt(end,5), plt2],{'Agent1','Agent2','Agent3','Agent4','Agent5','Centralized Filter'});
newPosition = [0.68 0.17 0.1 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',2);


%% Plotting the trace of the velocity components

MonteCarlo_TracePos(1).Data  = dlmread('Agent1_TraceVelCovariance_MonteCarlo.dat');
MonteCarlo_TracePos(2).Data  = dlmread('Agent2_TraceVelCovariance_MonteCarlo.dat');
MonteCarlo_TracePos(3).Data  = dlmread('Agent3_TraceVelCovariance_MonteCarlo.dat');
MonteCarlo_TracePos(4).Data  = dlmread('Agent4_TraceVelCovariance_MonteCarlo.dat');
MonteCarlo_TracePos(5).Data  = dlmread('Agent5_TraceVelCovariance_MonteCarlo.dat');

Centralized_TracePos(1).Data  = dlmread('Agent1_TraceVelCovariance_Centralized.dat');
Centralized_TracePos(2).Data  = dlmread('Agent2_TraceVelCovariance_Centralized.dat');
Centralized_TracePos(3).Data  = dlmread('Agent3_TraceVelCovariance_Centralized.dat');
Centralized_TracePos(4).Data  = dlmread('Agent4_TraceVelCovariance_Centralized.dat');
Centralized_TracePos(5).Data  = dlmread('Agent5_TraceVelCovariance_Centralized.dat');
figure('Renderer', 'painters', 'Position', [10 10 1000 700])
plt   = zeros(size(MonteCarlo_TracePos(1).Data,1), Num_deputies);
for k = 1:Num_deputies
    subplot(3,2,k);
    plt(:,k) = plot(indexY,MonteCarlo_TracePos(k).Data,'Color',ColorMatrix(k,:));
    hold on
    plt2     = plot(indexY,Centralized_TracePos(k).Data,'Color',c2);
    grid on
    xlabel('Period')
    ylabel(sprintf('Trace of P$_{%d}$', k))
    set(gca, 'YScale', 'log')
end
hL = legend([plt(end,1),plt(end,2),plt(end,3),plt(end,4),plt(end,5), plt2],{'Agent1','Agent2','Agent3','Agent4','Agent5','Centralized Filter'});
newPosition = [0.68 0.17 0.1 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',2);


%% Consistency Analysis
close all
clc
ConstAnalysis(1).Data  = dlmread('Agent1_ConstAnalysis_MonteCarlo.dat');
ConstAnalysis(2).Data  = dlmread('Agent2_ConstAnalysis_MonteCarlo.dat');
ConstAnalysis(3).Data  = dlmread('Agent3_ConstAnalysis_MonteCarlo.dat');
ConstAnalysis(4).Data  = dlmread('Agent4_ConstAnalysis_MonteCarlo.dat');
ConstAnalysis(5).Data  = dlmread('Agent5_ConstAnalysis_MonteCarlo.dat');
nn = 6; MCs = size(ConstAnalysis(5).Data,1)/nn;
Label = {'Agent1','Agent2','Agent3','Agent4','Agent5'};
figure('Renderer', 'painters', 'Position', [10 10 1000 700])
for ii = 1:MCs
    idx = 1+nn*(ii-1):nn*ii;
    for k = 1:Num_deputies
        plt   = zeros(6);
        subplot(n/2,2,k)
        plt2 = yline(3,'--','LineWidth', 2, 'Color', 'r');
        yline(-3,'--','LineWidth', 2, 'Color', 'r');
        hold on
        for i = 1:n
            plt(:,i) = plot(indexY,ConstAnalysis(k).Data(idx(i),:),'Color',ColorMatrix(i,:));
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


