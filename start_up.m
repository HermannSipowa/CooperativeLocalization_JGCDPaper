%% Figure setting

% to confirm the default settings
% get(groot,'factory')
%% font
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');


% GUI font
set(0, 'defaultUicontrolFontName', 'Meiryo UI');
% Axis font
set(0, 'defaultAxesFontName', 'Times');
% Title, legend font
set(0, 'defaultTextFontName', 'Times');

% GUI fontsize
set(0, 'defaultUicontrolFontSize', 9);

% Axos fontsize
set(0, 'defaultAxesFontSize', 18);
%set(0, 'defaultLabelFontSize', 15);
set(0, 'defaultLegendFontSize', 20);
%set(0, 'defaultPolarAxesFontSize', 18);
set(0, 'defaultTextarrowshapeFontSize', 18);
set(0, 'defaultTextboxshapeFontSize', 18);


% Title, legend fontsize
set(0, 'defaultTextFontSize', 20);

%% other
% axis line thickness
set(0, 'DefaultAxesLineWidth', 2);

% legend object line thickness
set(0, 'DefaultLineLineWidth', 2);

% default color map
cmap = spring(128);
set(0, 'defaultFigureColormap', cmap);

% plot color default
% corder = linspecer(9,'qualitative');
corder = [[0.3 1 0];[1 0 1];[1 0 0];[0 0 1];[0.7500 0 0.7500];[0.7500 0.7500 0];[0.2500 0.2500 0.2500];[1 0.7500 0.7500]];
set(0, 'defaultAxesColorOrder', corder);

clearvars cmap corder;