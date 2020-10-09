clear all
close all
clc
start_up
format short

color1 = rgb('DarkTurquoise');
color2 = rgb('DarkGreen');
Paa = [1.1287    1.0334;
    1.0334    1.3426];
Pbb = [1.6889    0.5746;
    0.5746    0.3636];

C1 = chol(Paa,'lower');
C2 = chol(Pbb,'lower');


fh = figure;
hold on
Determinant = 100;
N = 20;
omega = linspace(0.1,0.9,N);
rho   = linspace(0,0.9,N);
for i = 1:N
    
    Pab = rho(i)*C1*C2';
    Pba = rho(i)*C2*C1';
    Matrix = [Paa Pab;
        Pba Pbb];
    K1 = (Pbb - Pba)/(Paa + Pbb - Pab - Pba); K2 = eye(2) - K1;
    Pcc_known = [K1 K2] * Matrix * [K1 K2]';
    h3 = error_ellipse(Pcc_known);
    h3.Color = color1;
    h3.LineWidth = 1.5;
    
    
    
    k1 = omega(i); k2 = 1-omega(i);
    Pcc_Unknown = inv( k1*inv(Paa) + k2*inv(Pbb));
    h4 = error_ellipse(Pcc_Unknown);
    h4.Color = 'r';
    h4.LineWidth = 1;
    if trace(Pcc_Unknown)<Determinant
        Pcc_optimal = Pcc_Unknown;
        Determinant = trace(Pcc_Unknown);
    end
end

h1 = error_ellipse(Paa);
h1.Color = 'b';
h1.LineWidth = 1.5;

h2 = error_ellipse(Pbb);
h2.Color = color2;
h2.LineWidth = 1.5;

h5 = error_ellipse(Pcc_optimal);
h5.Color = 'k';
h5.Marker = 'none';
h5.LineWidth = 1.5;
% grid on
% ylim([-1.5 3])                 

% Hl = legend([h1, h2, h3],{'$P_{aa}$','$P_{bb}$','$P_{cc}$ (known $P_{ab}$)'},...
%     'Location','northwest','NumColumns',2);

Hl = legend([h1, h2, h3, h4, h5],{'$P_{aa}$','$P_{bb}$','$P_{cc}$ (known $P_{ab}$)',...
    '$P_{cc}$ (function of $\omega$)', '$P_{cc}$ (optimal $\omega^*$)'},'Location','northwest','NumColumns',2);
Hl.FontSize = 10;
newPosition = [0.47 0.81 0.1 0.1];
newUnits = 'normalized';
set(Hl,'Position', newPosition,'Units', newUnits,'NumColumns',2);

% set all units inside figure to normalized so that everything is scaling accordingly
set(findall(fh,'Units','pixels'),'Units','normalized');
% do not show figure on screen
% set(fh, 'visible', 'off')
% set figure units to pixels & adjust figure size
fh.Units = 'pixels';
fh.OuterPosition = [0 0 336 336];
% define resolution figure to be saved in dpi
res = 1000;
% recalculate figure size to be saved
set(fh,'PaperPositionMode','manual')
fh.PaperUnits = 'inches';
fh.PaperPosition = [0 0 4420 4420]/res;
% save figure
print(fh,'CovInter2D_UnKnownCorrelation','-dpng',sprintf('-r%d',res))

