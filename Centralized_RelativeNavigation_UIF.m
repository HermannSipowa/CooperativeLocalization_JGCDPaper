clear variables
close all
delete CentralizedData.mat
clc
start_up
format long e
mu_Earth = 3.986004415E5;
c1 = rgb('RosyBrown'); c2 = rgb('Black'); c3 = rgb('Lime');
c4 = rgb('Tomato'); c5 = rgb('DarkBlue'); c6 = rgb('DarkTurquoise');
ColorMatrix = [c1;c6;c3;c4;c5;c2];
dt = 60;

%% A) Defining the inital conditions for the spaceraft in the system
% Target initial conditions(Orbital Elements)
%*********************************************%
a1          = 1.5E+4; % Semi-major axis [Km] 7500; % 
e1          = 0.2;    % Eccentricity
inc1_deg    = 10;     % Inclination [deg]
BigOmg1_deg = 45;     % RAAN [deg]
LitOmg1_deg = 75;     % AOP [deg]
f_deg       = 0;      % True anomaly [deg]

% Converting angles in radians 
%******************************%
inc1 = deg2rad(inc1_deg); BigOmg1 = deg2rad(BigOmg1_deg);
LitOmg1 = deg2rad(LitOmg1_deg); f = deg2rad(f_deg);
COE1 = [a1,e1,inc1,BigOmg1,LitOmg1,f];
[Position_target,Velocity_target]  = COEstoRV(COE1,mu_Earth);
X_Chief0 = [Position_target; Velocity_target];
RelativeState = zeros(6,5);
for i = 1:5
    % Relative orbital elements of the deputy
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    dela      = 0;                  % Relative semi-major axis [Km]
    dele      = -2/(3*a1);          % Relative eccentricity
    deli      = i/(4*a1);           % Relative inclination [rad]
    delLitOmg = (-1) * i/(5*a1);    % Relative RAAN [rad]
    delBigOmg = (-1)^i * i/(5*a1);  % Relative AOP [rad]
    delf      = pi*i/(15*a1);       % Relative true anomaly [rad]
    
    
    % Chaser initial conditions(Orbital Elements)
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    a2          = a1 + dela;        % Semi-major axis [Km]
    e2          = e1 + dele;        % Eccentricity
    inc2_deg    = inc1_deg;         % Inclination  [deg]
    BigOmg2_deg = BigOmg1_deg;      % RAAN [rad]
    LitOmg2_deg = LitOmg1_deg;      % AOP [deg]
    f2_deg      = f_deg;            % True anomaly [deg]
    
    
    % Computing the ECI deputy position
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    inc2    = deg2rad(inc2_deg) + deli; BigOmg2 = deg2rad(BigOmg2_deg)+delBigOmg;
    LitOmg2 = deg2rad(LitOmg2_deg)+delLitOmg; f2 = deg2rad(f2_deg)+delf;
    COE2    = [a2,e2,inc2,BigOmg2,LitOmg2,f2];
    [Position_deputy,Velocity_deputy] = COEstoRV(COE2,mu_Earth);
    
    
    % Computing the relative states (as seen in the LVLH)
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    TN      = DCM(Position_target,Velocity_target); rt_norm = norm(Position_target);
    h_vec   = cross(Position_target,Velocity_target); h_norm = norm(h_vec);
    eh      = h_vec/h_norm;
    N_nudot = h_vec/rt_norm^2;
    NR_rel  = Position_deputy - Position_target; NV_rel2 = Velocity_deputy - Velocity_target;
    TR_rel0 = TN*NR_rel; TV_rel0 = TN*(NV_rel2 - cross(N_nudot,NR_rel));
    
    % Storing the relative states
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    RelativeState(:,i) = [TR_rel0; TV_rel0];
    
end

%% B) Computing all the true trajectories
% Setting the integrator parameters
%***********************************%
[n, Num_deputies] = size(RelativeState);
Period            = 2*pi*sqrt(a1^3/mu_Earth);
IntegrationTime   = 6*Period;
tvec              = 0:dt:IntegrationTime;
options           = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-20);

% Intagrading the chief's and deputies' trajectory
%**************************************************%
RelativeState = [X_Chief0, RelativeState];
[~, Xrel]     = ode113(@(t,X_aug)UnperturbedRelativeMotionODE(t,X_aug,mu_Earth,n),tvec,RelativeState,options);
X_Chief       = Xrel(:,1:n);
Xrel          = Xrel(:,n+1:end);
index         = 1:6;
X_Deputy      = struct('x', cell(1, Num_deputies));
for k = 1:Num_deputies
    X_Deputy(k).x = Xrel(:,index)';
    index         = index + 6;
end

%% C) Visualazing the geomety of the problem
% Trajectory in the chief's LVLH frame
%**************************************%
for j = 1
% figure
% Label = {'Agent1','Agent2','Agent3','Agent4','Agent5'};
% plt   = zeros(5);
% plot3(0,0,0,'bo',...
%     'LineWidth',2,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor',c2,...
%     'MarkerSize',15);
% hold on
% grid on
% xlabel('X[km]')
% ylabel('Y[km]')
% zlabel('Z[km]')
% 
% for k = 1:Num_deputies
%     plt(k) = plot3(X_Deputy(k).x(1,:),X_Deputy(k).x(2,:),X_Deputy(k).x(3,:));
%     hold on
%     plot3(X_Deputy(k).x(1,1),X_Deputy(k).x(2,1),X_Deputy(k).x(3,1),'b*',...
%     'LineWidth',2,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor',c1,...
%     'MarkerSize',5);
%     grid on
%     legend(plt(1:k),Label(1:k))
% end
% axis equal
end

%% D) Generating the noisy measurements
p = 4;
Yrel = nan(p+2, Num_deputies, length(tvec));
Yabsolute = nan(p+3, length(tvec));
sigma1_measurement = 1E-3;
sigma2_measurement = 1E-6;
sigma3_measurement = 1E-3;
sigma4_measurement = 1E-3;

% Absolute measurements from the chief and ground sensors
%*********************************************************%
origin = zeros(6,1);
sensor3 = [1000*ones(3,1); zeros(3,1)];
T_delay = dt; count = 0;
for i = 1:length(tvec)
%     if floor(tvec(i)/T_delay)>=count
%         index = randi([4 5]);
%         sensorID = randi([0 2]);
%         coef = randn(p,1);
%         v = [sigma1_measurement; sigma2_measurement; sigma3_measurement; sigma4_measurement].*coef;
%         Xi = X_Deputy(index).x(:,i);
%         if sensorID == 0
%             Xj = origin;
%         elseif sensorID == 1
%             Xj = X_Chief(i,:)';
%         elseif sensorID == 2
%             Xj = sensor3;
%         end
%         count = count+1;
%         Yabsolute(:, i) = [tvec(i); MeasurementFunc(Xi,Xj) + v; index; sensorID];
%     end
end

% Relative measurements between deputies
%****************************************%
for i = 1:length(tvec)
    target = randi([1 Num_deputies]);
    % index  = randperm(Num_deputies, Num_deputies-1);
    index  = randperm(Num_deputies,randi(Num_deputies-1));
    
    if numel(find(index==target)) ~= 0
        while length(index) ~= length(unique(index)) || numel(find(index==target)) ~= 0
            index(index==target) = randi(Num_deputies);
            [~, w] = unique( index, 'stable' );
            ind = setdiff( 1:numel(index), w );
            if ~isempty(ind)
                index(ind) = randi(Num_deputies);
            end
        end
    end
    
    for k = 1:length(index)
        coef = randn(p,1);
        v = [sigma1_measurement; sigma2_measurement; sigma3_measurement; sigma4_measurement].*coef;
        Xi = X_Deputy(target).x(:,i);
        Xj = X_Deputy(index(k)).x(:,i);
        Yrel(:, index(k), i) = [tvec(i); MeasurementFunc(Xi,Xj) + v; target];
    end
end


%% (1) Initializing the UKF
Cov = diag([(1E-1)^2, (1E-1)^2, (1E-1)^2, (1E-5)^2, (1E-5)^2, (1E-5)^2]);
Q0  = diag([(1E-5)^2, (1E-5)^2, (1E-5)^2, (1E-6)^2, (1E-6)^2, (1E-6)^2]);
for k = 1:Num_deputies
    if k == 1
        P_post = Cov;
        X_post = X_Deputy(k).x(:,1) + Cov*randn(n,1);
        Q = Q0;
    else
        P_post = blkdiag(P_post,Cov);
        X_post = [X_post; X_Deputy(k).x(:,1) + Cov*randn(n,1);];
        Q = blkdiag(Q,Q0);
    end
end
Q = zeros(n*Num_deputies);
R = diag([(sigma1_measurement)^2 (sigma2_measurement)^2 ...
    sigma3_measurement^2 sigma4_measurement^2]); % Measurement covariance


%% (2) Calculating all the weights
alpha = 1; % Alpha varies from 1E-4 to 1
beta  = 2;
L     = n*Num_deputies; % Number of state to be estimated
kappa = 3-L;
lamda = alpha^2*(L+kappa)-L;
gamma = sqrt(L+lamda);

Wo_m = lamda/(L+lamda);
Wo_c = lamda/(L+lamda)+(1-alpha^2+beta);
Wi = zeros(1,2*L);
for i = 1:2*L
    Wi(i) = 1/(2*(L+lamda));
end
Wm = [Wo_m, Wi];
Wc = [Wo_c, Wi];


m = length(tvec);
progressbar('Simulation') % Init 1 progress bar
for epoch = 1:m % Measurement start at t=0
    
    ti = tvec(epoch);
    progressbar(epoch/m) % Update the progress bar
    
    %% (3) Propagating the state of all the agents at once
    % (3.a) Computing Sigma Points (matrix) for t-1
    %**********************************************%
    S_post = sqrtm(P_post);
    Chi_post = [X_post, X_post + gamma*S_post, X_post - gamma*S_post];
    Num_trajs = size(Chi_post,2);
    
    % (3.b) Propagate Sigma Points through non-linear system
    %*******************************************************%
    if ti~=0
        for kk = 1:Num_trajs
            if kk == 1
                Chi_post_reshaped =  reshape(Chi_post(:,kk),6,5);
            else
                Chi_post_reshaped = [Chi_post_reshaped, reshape(Chi_post(:,kk),6,5)];
            end
        end
        X_aug0 = [X_Chief(epoch-1,:)', Chi_post_reshaped]; [n,mm] = size(X_aug0); index = 1:5;
        [~, Xtrajs] = ode113(@(t,X_aug)UnperturbedRelativeMotionODE(t,X_aug,mu_Earth,n),[0 dt],X_aug0,options);
        Xtrajs = reshape(Xtrajs(end,:),[n,mm]); Xtrajs = Xtrajs(:,2:end);
        for kk = 1:Num_trajs
            if kk == 1
                Xnom = reshape(Xtrajs(:,index),30,1);
            else
                Xnom = [Xnom, reshape(Xtrajs(:,index),30,1)];
            end
            index = index + 5;
        end
    else
        Xnom = Chi_post;
    end
    
    % (3.c) Perform a state and covariance time update
    %**************************************************%
    X_ap = zeros(L,1);
    P_ap = zeros(L,L);
    for j = 1:Num_trajs
        Xnomi = Xnom(:,j);
        X_ap = X_ap + Wm(j)*Xnomi; % mean of weighted sigma points
    end
    
    for j = 1:Num_trajs
        Xnomi = Xnom(:,j);
        P_ap = P_ap + Wc(j)*(Xnomi-X_ap)*(Xnomi-X_ap)'; % covariance of weighted sigma points
    end
    P_ap = Q + P_ap;
    
    % (3.d) Re-compute Sigma Points to incorporate effects of process noise
    %***********************************************************************%
    S_ap = sqrtm(P_ap);
    Chi_ap = [X_ap, X_ap + gamma*S_ap, X_ap - gamma*S_ap];
    clear P_post X_post

    

    %% (4) Message Creation and Sending
    Y_post = zeros(L,L); z_post = zeros(L,1);
    for jj = 1:Num_deputies % The index "jj" refers to the agent doing the tracking
        idx_Sensor = 1+n*(jj-1):n*jj;
        for ii = 1:Num_deputies % The index "ii" refers to the agent being tracked
            idx_Agent = 1+n*(ii-1):n*ii;
            y_target = zeros(p,Num_trajs); ys_sigma = zeros(p,Num_trajs);
            y = zeros(p,1); y_ap = zeros(p,1); ys_ap = zeros(p,1);
            if Yrel(end, jj, epoch) == ii % Checking if the sensor "jj" measured of agent "ii"
                y = Yrel(2:p+1, jj, epoch);
                for j = 1:Num_trajs
                    Xi = Chi_ap(idx_Agent,j);
                    Xj = X_Deputy(jj).x(:,epoch); % X_ap(idx_Sensor); % 
                    y_target(:,j) = MeasurementFunc(Xi,Xj);
                    
                    Xii = X_ap(idx_Agent);
                    Xjj = Chi_ap(idx_Sensor,j);
                    ys_sigma(:,j) = MeasurementFunc(Xii,Xjj);
                end
            end
            
            for j = 1:Num_trajs
                yi    = y_target(:,j);
                y_ap  = y_ap+Wm(j)*yi;
                
                yii   = ys_sigma(:,j);
                ys_ap = ys_ap+Wm(j)*yii;
            end
            
            
            % (4.b) Compute the innovation and cross-correlation covariances
            %***************************************************************%
            nn    = size(y_ap,1);
            P_xy  = zeros(L,nn);
            Ps_xy = zeros(L,nn);
            
            for j = 1:Num_trajs
                Xnomi  = Chi_ap(:,j);
                yi     = y_target(:,j); yii = ys_sigma(:,j);
                P_xy   = P_xy  + Wc(j)*(Xnomi-X_ap)*(yi-y_ap)';
                Ps_xy  = Ps_xy + Wc(j)*(Xnomi-X_ap)*(yii-ys_ap)';
            end
            
            
            % (4.c) Message sending/reception
            %*********************************%
            H = (P_ap\P_xy)'; Hs = (P_ap\Ps_xy)';
            v = y - y_ap;
            % if det(Hs*P_ap*Hs') ~= 0
            %     det(Hs*P_ap*Hs')
            % end
            
            z_post = z_post + H'/(R + Hs*P_ap*Hs')*(v + H*X_ap); % + Hs*P_ap*Hs'
            Y_post = Y_post + H'/(R + Hs*P_ap*Hs')*H;
        end
    end
    
    %% (5) Umpdating the statellite state
    % (5.1) Computing the information vector and matrix
    %***************************************************%
    z_post = P_ap\X_ap + z_post; % information vector
    Y_post = inv(P_ap) + Y_post; % information matrix
    
    % (5.1) Computing the posteriori state and covarinace
    %*****************************************************%
    X_post = Y_post\z_post;
    P_post = inv(Y_post);
    
    clear P_ap X_ap Chi_ap Y_post z_post
   
    %% (6) Computing the measurement residuals
    % Computing the residuals for relative measurements
    %**************************************************%
    for k = 1:Num_deputies
        % Storing estimated state
        %*************************%
        idx_Agent = 1+n*(k-1):n*k;
        X_updated(k).x(:,epoch) = X_post(idx_Agent);
        
        % Matrix allocation
        %*******************%
        PostfitAgent1(k).y(:,epoch) = nan(p,1);
        PostfitAgent2(k).y(:,epoch) = nan(p,1);
        PostfitAgent3(k).y(:,epoch) = nan(p,1);
        PostfitAgent4(k).y(:,epoch) = nan(p,1);
        PostfitAgent5(k).y(:,epoch) = nan(p,1);
        
        % Computing the residuals
        %*************************%
        for j = 1:Num_deputies % Referening to the agent doing the tracking
            if j==k
                continue
            end
            
            if Yrel(end,j,epoch)==k % Check if agent j took a measurement about agent k
                idx_Sensor  = 1+n*(j-1):n*j;
                Xi = X_post(idx_Agent);
                Xj = X_Deputy(j).x(:,epoch);
                G = MeasurementFunc(Xi,Xj);
                y = Yrel(2:p+1,j,epoch);
                % y - G
                if j==1
                    PostfitAgent1(k).y(:,epoch) = y - G;
                elseif j==2
                    PostfitAgent2(k).y(:,epoch) = y - G;
                elseif j==3
                    PostfitAgent3(k).y(:,epoch) = y - G;
                elseif j==4
                    PostfitAgent4(k).y(:,epoch) = y - G;
                elseif j==5
                    PostfitAgent5(k).y(:,epoch) = y - G;
                end
            end
        end
        % Sampeling the covariance inthe error
        %**************************************%
        Sigma(k).s(:,epoch)  = 3*sqrt(diag(P_post(idx_Agent,idx_Agent)));
        TracePos(k).t(epoch) = trace(P_post(idx_Agent(1:3),idx_Agent(1:3)));
        TraceVel(k).t(epoch) = trace(P_post(idx_Agent(4:6),idx_Agent(4:6)));
    end
    
end
% Storing the trace of the convariance matrice for every element
PosData = [TracePos(1).t(:)';TracePos(2).t(:)';TracePos(3).t(:)';...
      TracePos(4).t(:)';TracePos(5).t(:)'];
save('CentralizedData.mat','PosData')
Mat = matfile('CentralizedData.mat','Writable',true);
VelData = [TraceVel(1).t(:)';TraceVel(2).t(:)';TraceVel(3).t(:)';...
      TraceVel(4).t(:)';TraceVel(5).t(:)'];
Mat.VelData = VelData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%
% Plotting the measurements residuals
Hist_index   = 1:m; %find(tvec >= Period,1);
indexY       = tvec(1:m)/Period;
sigma_rho    = 3*sqrt(R(1,1))*ones(1,m);  % Sampeling the covariance in the error
sigma_rhodot = 3*sqrt(R(2,2))*ones(1,m);  % Sampeling the covariance in the error
sigma_El     = 3*sqrt(R(3,3))*ones(1,m);  % Sampeling the covariance in the error
sigma_Az     = 3*sqrt(R(4,4))*ones(1,m);  % Sampeling the covariance in the error
Label2       = {'Range[km]', 'Range-Rate[km/s]', 'El[rad]', 'Az[rad]'};
Label3       = {'Range Error','Range-Rate Error', 'El Error', 'Az Error'};
Sigma_rho    = [sigma_rho;sigma_rhodot;sigma_El;sigma_Az];


for k = 1%:Num_deputies
    figure('Renderer', 'painters', 'Position', [10 10 900 600])
    for i=1:4
        if i == 2
            iii = 3;
        elseif i == 3
            iii = 5;
        elseif i == 4
            iii = 7;
        else
            iii = i;
        end
        subplot(4,2,iii)
        plt1 = plot(indexY,PostfitAgent1(k).y(i,:),'o','Color', c1);
        hold on
        plt2 = plot(indexY,PostfitAgent2(k).y(i,:),'o','Color', c2);
        plt3 = plot(indexY,PostfitAgent3(k).y(i,:),'o','Color', c3);
        plt4 = plot(indexY,PostfitAgent4(k).y(i,:),'o','Color', c4);
        plt5 = plot(indexY,PostfitAgent5(k).y(i,:),'o','Color', c5);
        plot(indexY,Sigma_rho(i,:),'r--',indexY,-Sigma_rho(i,:),'r--')
        hold off
        ylabel(Label2(i))
        xlabel('Period')
        % legend([plt1,plt2,plt3,plt4,plt5],{'Agent1','Agent2','Agent3','Agent4','Agent5'})
        
        
        subplot(4,2,iii+1)
        plt1 = histogram(PostfitAgent1(k).y(i,Hist_index:end),'FaceColor',c1);
        hold on
        plt2 = histogram(PostfitAgent2(k).y(i,Hist_index:end),'FaceColor',c2);
        plt3 = histogram(PostfitAgent3(k).y(i,Hist_index:end),'FaceColor',c3);
        plt4 = histogram(PostfitAgent4(k).y(i,Hist_index:end),'FaceColor',c4);
        plt5 = histogram(PostfitAgent5(k).y(i,Hist_index:end),'FaceColor',c5);
        xline(Sigma_rho(i,1),'--','LineWidth', 2, 'Color', 'r');
        xline(-Sigma_rho(i,1),'--','LineWidth', 2, 'Color', 'r');
        xlabel(Label3(i))
        set(gca,'view',[90 -90])
    end
    % sgt = sgtitle(['Relative Measurements Residual (Agent' num2str(k) ')']);
    % sgt.FontSize = 35;
    
end

%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Plotting the StateError Residuals
Label={'$x$[km]', '$\dot{x}$[m/s]', '$y$[km]', '$\dot{y}$[m/s]', '$z$[km]', '$\dot{z}$[m/s]'};
for k = 1%:Num_deputies
    StateError = X_updated(k).x - X_Deputy(k).x(:,1:m);
    StateError = [StateError(1,:);1E3*StateError(4,:);StateError(2,:);...
                  1E3*StateError(5,:);StateError(3,:);1E3*StateError(6,:)];
    Sigma(k).s = [Sigma(k).s(1,:);1E3*Sigma(k).s(4,:);Sigma(k).s(2,:);...
                  1E3*Sigma(k).s(5,:);Sigma(k).s(3,:);1E3*Sigma(k).s(6,:)];
    figure('Renderer', 'painters', 'Position', [10 10 900 600])
    for i=1:n
        subplot(n/2,2,i)
        plot(indexY,StateError(i,1:m),'b')
        hold on
        plot(indexY,Sigma(k).s(i,:),'r--',indexY,-Sigma(k).s(i,:),'r--')
        hold off
        ylabel(Label(i))
        xlabel('Period')
    end
    % Construct a Legend with the data from the sub-plots
    hL = legend('State Error$~~~~~~~~~~~~~~~$','3$\sigma$ Covariance Envelope');
    % Programatically move the Legend
    newPosition = [.27 .94 0.5 0.05];
    newUnits = 'normalized';
    set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',2);
    % sgt = sgtitle(['State Residual (Agent' num2str(k) ')']);
    % sgt.FontSize = 35;
    
end

%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Plotting the trace of the covariance matrices
figure
Label = {'Agent1','Agent2','Agent3','Agent4','Agent5'};
plt   = zeros(5);
for k = 1:Num_deputies
    plt(k) = plot(indexY,TracePos(k).t(:));
    hold on
    grid on
    legend(plt(1:k),Label(1:k))
end
xlabel('Period')
ylabel('Trace of P')
set(gca, 'YScale', 'log')

% %%
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% % Plotting the trace of the covariance matrices
% figure
% Label = {'Agent1','Agent2','Agent3','Agent4','Agent5'};
% plt   = zeros(5);
% for k = 1:Num_deputies
%     plt(k) = plot(indexY,TraceVel(k).t(:));
%     hold on
%     grid on
%     legend(plt(1:k),Label(1:k))
% end
% xlabel('Period')
% ylabel('Trace of P')
% set(gca, 'YScale', 'log')

