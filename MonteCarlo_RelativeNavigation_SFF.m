clear all
close all
% delete Agent1_TracePosCovariance_MonteCarlo.dat
% delete Agent2_TracePosCovariance_MonteCarlo.dat
% delete Agent3_TracePosCovariance_MonteCarlo.dat
% delete Agent4_TracePosCovariance_MonteCarlo.dat
% delete Agent5_TracePosCovariance_MonteCarlo.dat
% delete Agent1_TraceVelCovariance_MonteCarlo.dat
% delete Agent2_TraceVelCovariance_MonteCarlo.dat
% delete Agent3_TraceVelCovariance_MonteCarlo.dat
% delete Agent4_TraceVelCovariance_MonteCarlo.dat
% delete Agent5_TraceVelCovariance_MonteCarlo.dat
% delete Agent1_ConstAnalysis_MonteCarlo.dat
% delete Agent2_ConstAnalysis_MonteCarlo.dat
% delete Agent3_ConstAnalysis_MonteCarlo.dat
% delete Agent4_ConstAnalysis_MonteCarlo.dat
% delete Agent5_ConstAnalysis_MonteCarlo.dat
clc
system('caffeinate -dims &');
start_up
format long e
mu_Earth = 3.986004415E5;
c1 = rgb('RosyBrown'); c2 = rgb('Black'); c3 = rgb('Lime');
c4 = rgb('Tomato'); c5 = rgb('DarkBlue'); c6 = rgb('DarkTurquoise');
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
RelativeState = [];
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
    a2          = a1 + dela;       % Semi-major axis [Km]
    e2          = e1 + dele;       % Eccentricity
    inc2_deg    = inc1_deg;        % Inclination  [deg]
    BigOmg2_deg = BigOmg1_deg;     % RAAN [rad]
    LitOmg2_deg = LitOmg1_deg;     % AOP [deg]
    f2_deg      = f_deg;           % True anomaly [deg]
    
    
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
    RelativeState = [RelativeState, [TR_rel0; TV_rel0]];
    
end

%% B) Computing all the true trajectories
% Setting the integrator parameters
%***********************************%
[n, Num_deputies] = size(RelativeState);
Period            = 2*pi*sqrt(a1^3/mu_Earth);
IntegrationTime   = 3*Period;
tvec              = 0:dt:IntegrationTime;
options           = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-20);

% Intagrading the chief's and deputies' trajectory
%**************************************************%
RelativeState = [X_Chief0, RelativeState];
[~, Xrel]     = ode113(@(t,X_aug)UnperturbedRelativeMotionODE(t,X_aug,mu_Earth,n),tvec,RelativeState,options);
X_Chief       = Xrel(:,1:n);
Xrel          = Xrel(:,n+1:end);
index         = 1:6;
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



%% Monte Carlo
MC_runs = 20;
progressbar('Monte Carlo Trials','Simulation') % Init 2 bars
for monte = 1:MC_runs
    
   progressbar([],0) % Reset 2nd bar
    
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
    origin  = zeros(6,1);
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
    for k = 1:Num_deputies
        P_post(k).P = Cov;
        X_post(k).x = X_Deputy(k).x(:,1) + Cov*randn(n,1);
    end
    Q = zeros(n);
    R = diag([(sigma1_measurement)^2 (sigma2_measurement)^2 ...
        sigma3_measurement^2 sigma4_measurement^2]); % Measurement covariance
    
    
    %% (2) Calculating all the weights
    alpha = 1; % Alpha varies from 1E-4 to 1
    beta  = 2;
    L     = n; % Number of state to be estimated
    kappa = 3-L;
    lamda = alpha^2*(L+kappa)-L;
    gamma = sqrt(L+lamda);
    
    Wo_m = lamda/(L+lamda);
    Wo_c = lamda/(L+lamda)+(1-alpha^2+beta);
    Wi = zeros(1,12);
    for i = 1:2*L
        Wi(i) = 1/(2*(L+lamda));
    end
    Wm = [Wo_m, Wi];
    Wc = [Wo_c, Wi];
    
   m = length(tvec);
for epoch = 1:m % Measurement start at t=0
    
    ti = tvec(epoch);
    progressbar([],epoch/m) % Update 2rd bar
    
    %% (3) Separately propagate the agents states
    for k = 1:Num_deputies
        % (3.a) Computing Sigma Points (matrix) for t-1
        %**********************************************%
        S_post = sqrtm(P_post(k).P);
        Chi_post = [X_post(k).x, X_post(k).x + gamma*S_post, X_post(k).x - gamma*S_post];
        Num_worker = size(Chi_post,2);
        
        % (3.b) Propagate Sigma Points through non-linear system
        %*******************************************************%
        if ti~=0
            X_aug0    = [X_Chief(epoch-1,:)', Chi_post]; [n,mm] = size(X_aug0);
            [~, Xnom] = ode113(@(t,X_aug)UnperturbedRelativeMotionODE(t,X_aug,mu_Earth,n),[0 dt],X_aug0,options);
            Xnom      = reshape(Xnom(end,:),[n,mm]); Xnom = Xnom(:,2:end);
        else
            Xnom      = Chi_post;
        end
        
        % (3.c) Perform a state and covariance time update
        %**************************************************%
        StuctX_ap(k).x  = zeros(n,1);
        StructP_ap(k).P = zeros(n,n);
        for j = 1:Num_worker
            Xnomi = Xnom(:,j);
            StuctX_ap(k).x = StuctX_ap(k).x + Wm(j)*Xnomi; % mean of weighted sigma points
        end
        
        for j = 1:Num_worker
            Xnomi = Xnom(:,j);
            StructP_ap(k).P = StructP_ap(k).P + Wc(j)*...
                (Xnomi-StuctX_ap(k).x)*(Xnomi-StuctX_ap(k).x)'; % covariance of weighted sigma points
        end
        StructP_ap(k).P = Q + StructP_ap(k).P;
        
        % (3.d) Re-compute Sigma Points to incorporate effects of process noise
        %***********************************************************************%
        S_ap = sqrtm(StructP_ap(k).P);
        StructChi_ap(k).x = [StuctX_ap(k).x,...
                             StuctX_ap(k).x + gamma*S_ap,...
                             StuctX_ap(k).x - gamma*S_ap];  
    end
    clear P_post X_post


    %% (4) Message Creation and Sending
    for jj = 1:Num_deputies % Index "jj" is referening to the agent doing the tracking
        Chi_ap_sensor = StructChi_ap(jj).x;  % StructChi_dist(jj).x; 
        X_sensor      = StuctX_ap(jj).x;     % StuctX_dist(jj).x;
        P_ap_sensor   = StructP_ap(jj).P;    % StructP_dist(jj).P;
        
        for ii = 1:Num_deputies % index "ii" is referening the agent being tracked
            Chi_ap   = StructChi_ap(ii).x; % StructChi_dist(ii).x;
            X_ap     = StuctX_ap(ii).x;    % StuctX_dist(ii).x;
            P_ap     = StructP_ap(ii).P;   % StructP_dist(ii).P;
            y_target = zeros(p,Num_worker); y = zeros(p,1);
            ys_sigma = zeros(p,Num_worker); ys_ap = zeros(p,1);
            if Yrel(end, jj, epoch) == ii % Checking if agent "jj" took a measurement of agent "ii"
                y = Yrel(2:p+1, jj, epoch);
                for j = 1:Num_worker
                    Xi = Chi_ap(:,j);
                    Xj = X_sensor;
                    y_target(:,j) = MeasurementFunc(Xi,Xj);
                    
                    Xii = X_ap;
                    Xjj = Chi_ap_sensor(:,j);
                    ys_sigma(:,j) = MeasurementFunc(Xii,Xjj);
                end
            end
            
            y_ap = zeros(p,1);
            for j = 1:Num_worker
                yi    = y_target(:,j);
                y_ap  = y_ap+Wm(j)*yi;
                
                yii   = ys_sigma(:,j);
                ys_ap = ys_ap+Wm(j)*yii;
            end
            
            % (4.b) Compute the innovation and cross-correlation covariances
            %***************************************************************%
            nn = size(y_ap,1);
            P_xy = zeros(n,nn);
            Ps_xy = zeros(n,nn);
            
            for j = 1:Num_worker
                Xnomi  = Chi_ap(:,j);
                yi     = y_target(:,j);
                P_xy   = P_xy + Wc(j)*(Xnomi-X_ap)*(yi-y_ap)';
                
                Xnomii = Chi_ap_sensor(:,j);
                yii    = ys_sigma(:,j);
                Ps_xy  = Ps_xy + Wc(j)*(Xnomii-X_sensor)*(yii-ys_ap)';
            end
            
            % (4.c) Message sending/reception
            %*********************************%
            H = (P_ap\P_xy)'; Hs = (P_ap_sensor\Ps_xy)';
            v = y - y_ap;
            
            FromNeighbors(jj,ii).ivector = H'/(R + Hs*P_ap_sensor*Hs')*(v + H*X_ap); 
            FromNeighbors(jj,ii).Imatrix = H'/(R + Hs*P_ap_sensor*Hs')*H;            
            
        end
    end
    %% (5) Message Absorption via Covariance Intersection
    % (5.1) Performing covatiance intersection on all the estimates
    %***************************************************************%
    for ii = 1:Num_deputies
        CIFused_Info.ivector = StructP_ap(ii).P\StuctX_ap(ii).x;
        CIFused_Info.Imatrix = inv(StructP_ap(ii).P);
        for jj = 1:Num_deputies
            currInfoMessage_ji   = FromNeighbors(jj,ii);
            [CIFused_Info,omega,fmin] = CI(CIFused_Info,currInfoMessage_ji);
        end
        
        X_post(ii).x = CIFused_Info.Imatrix\CIFused_Info.ivector; 
        P_post(ii).P = inv(CIFused_Info.Imatrix); 
    end
    

    
    %% (7) Computing the measurement residuals
    % Matrix allocation
    %*******************%
    AbsPostfitSensor0.y(:,epoch) = nan(p,1);
    AbsPostfitSensor1.y(:,epoch) = nan(p,1);
    AbsPostfitSensor2.y(:,epoch) = nan(p,1);
    % Computing the residuals for absolute measurements
    %**************************************************%
    for k = 5 % Referening to the agent doing the tracking
        if Yabsolute(end-1, epoch)==k % Check if agent k took a measurement about agent kk
            if Yabsolute(end,epoch) == 0
                Xj = origin;
            elseif Yabsolute(end,epoch) == 1
                Xj = X_Chief(epoch,:)';
            elseif Yabsolute(end,epoch) == 2
                Xj = sensor3;
            end
            Xi = X_post(k).x;
            G = MeasurementFunc(Xi,Xj);
            y = Yabsolute(2:p+1,epoch);
            
            if Yabsolute(end,epoch) == 0
                AbsPostfitSensor0.y(:,epoch) = y - G;
            elseif Yabsolute(end,epoch) == 1
                AbsPostfitSensor1.y(:,epoch) = y - G;
            elseif Yabsolute(end,epoch) == 2
                AbsPostfitSensor2.y(:,epoch) = y - G;
            end
        end
    end
    
    for k = 1:Num_deputies
        
        % Storing estimated state
        %*************************%
        X_updated(k).x(:,epoch) = X_post(k).x;
        
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
            if Yrel(end,j,epoch)==k % Check if agent k took a measurement about agent kk
                Xi = X_post(k).x;
                Xj = X_Deputy(j).x(:,epoch); % X_post(j).x; % 
                G = MeasurementFunc(Xi,Xj);
                y = Yrel(2:p+1,j,epoch);
                
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
        Sigma(k).s(:,epoch) = 3*sqrt(diag(P_post(k).P));
        % Trace(k).t(epoch)   = det(P_post(k).P);
        TracePos(k).t(epoch) = trace(P_post(k).P(1:3,1:3));
        TraceVel(k).t(epoch) = trace(P_post(k).P(4:6,4:6));
    end
    
    clear P_ap X_ap Chi_ap CIFused_Info FromNeighbors
end

% Storing data for the consistency monte carlo simulation
    for k = 1:Num_deputies
        StateError = X_updated(k).x - X_Deputy(k).x(:,1:m);
        StateError = [StateError(1,:); StateError(2,:); StateError(3,:); StateError(4,:); StateError(5,:); StateError(6,:)];
        ConstAnalysis(k).Data = StateError./Sigma(k).s;
    end
    % Storing the trace of the convariance matrice for every element
    if monte == 1
        dlmwrite('Agent1_ConstAnalysis_MonteCarlo.dat',ConstAnalysis(1).Data)
        dlmwrite('Agent2_ConstAnalysis_MonteCarlo.dat',ConstAnalysis(2).Data)
        dlmwrite('Agent3_ConstAnalysis_MonteCarlo.dat',ConstAnalysis(3).Data)
        dlmwrite('Agent4_ConstAnalysis_MonteCarlo.dat',ConstAnalysis(4).Data)
        dlmwrite('Agent5_ConstAnalysis_MonteCarlo.dat',ConstAnalysis(5).Data)
        
        % dlmwrite('Agent1_TracePosCovariance_MonteCarlo.dat',TracePos(1).t(:)')
        % dlmwrite('Agent2_TracePosCovariance_MonteCarlo.dat',TracePos(2).t(:)')
        % dlmwrite('Agent3_TracePosCovariance_MonteCarlo.dat',TracePos(3).t(:)')
        % dlmwrite('Agent4_TracePosCovariance_MonteCarlo.dat',TracePos(4).t(:)')
        % dlmwrite('Agent5_TracePosCovariance_MonteCarlo.dat',TracePos(5).t(:)')
        % dlmwrite('Agent1_TraceVelCovariance_MonteCarlo.dat',TraceVel(1).t(:)')
        % dlmwrite('Agent2_TraceVelCovariance_MonteCarlo.dat',TraceVel(2).t(:)')
        % dlmwrite('Agent3_TraceVelCovariance_MonteCarlo.dat',TraceVel(3).t(:)')
        % dlmwrite('Agent4_TraceVelCovariance_MonteCarlo.dat',TraceVel(4).t(:)')
        % dlmwrite('Agent5_TraceVelCovariance_MonteCarlo.dat',TraceVel(5).t(:)')
    else
        dlmwrite('Agent1_ConstAnalysis_MonteCarlo.dat',ConstAnalysis(1).Data,'-append')
        dlmwrite('Agent2_ConstAnalysis_MonteCarlo.dat',ConstAnalysis(2).Data,'-append')
        dlmwrite('Agent3_ConstAnalysis_MonteCarlo.dat',ConstAnalysis(3).Data,'-append')
        dlmwrite('Agent4_ConstAnalysis_MonteCarlo.dat',ConstAnalysis(4).Data,'-append')
        dlmwrite('Agent5_ConstAnalysis_MonteCarlo.dat',ConstAnalysis(5).Data,'-append')
        
        % dlmwrite('Agent1_TracePosCovariance_MonteCarlo.dat',TracePos(1).t(:)','-append')
        % dlmwrite('Agent2_TracePosCovariance_MonteCarlo.dat',TracePos(2).t(:)','-append')
        % dlmwrite('Agent3_TracePosCovariance_MonteCarlo.dat',TracePos(3).t(:)','-append')
        % dlmwrite('Agent4_TracePosCovariance_MonteCarlo.dat',TracePos(4).t(:)','-append')
        % dlmwrite('Agent5_TracePosCovariance_MonteCarlo.dat',TracePos(5).t(:)','-append')
        % dlmwrite('Agent1_TraceVelCovariance_MonteCarlo.dat',TraceVel(1).t(:)','-append')
        % dlmwrite('Agent2_TraceVelCovariance_MonteCarlo.dat',TraceVel(2).t(:)','-append')
        % dlmwrite('Agent3_TraceVelCovariance_MonteCarlo.dat',TraceVel(3).t(:)','-append')
        % dlmwrite('Agent4_TraceVelCovariance_MonteCarlo.dat',TraceVel(4).t(:)','-append')
        % dlmwrite('Agent5_TraceVelCovariance_MonteCarlo.dat',TraceVel(5).t(:)','-append')
    end
    progressbar(monte/MC_runs) % Update 1st bar
end % Ending the Monte Carlo loop
system('killall caffeinate');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%
% Plotting the measurements residuals
Hist_index   = find(tvec >= Period,1);
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
        line1 = plot(indexY,PostfitAgent1(k).y(i,:),'o','Color', c1);
        hold on
        line2 = plot(indexY,PostfitAgent2(k).y(i,:),'o','Color', c2);
        line3 = plot(indexY,PostfitAgent3(k).y(i,:),'o','Color', c3);
        line5 = plot(indexY,PostfitAgent4(k).y(i,:),'o','Color', c4);
        plt5 = plot(indexY,PostfitAgent5(k).y(i,:),'o','Color', c5);
        plot(indexY,Sigma_rho(i,:),'r--',indexY,-Sigma_rho(i,:),'r--')
        hold off
        ylabel(Label2(i))
        xlabel('Period')
        % legend([plt1,plt2,plt3,plt4,plt5],{'Agent1','Agent2','Agent3','Agent4','Agent5'})
        
        
        subplot(4,2,iii+1)
        line1 = histogram(PostfitAgent1(k).y(i,Hist_index:end),'FaceColor',c1);
        hold on
        line2 = histogram(PostfitAgent2(k).y(i,Hist_index:end),'FaceColor',c2);
        line3 = histogram(PostfitAgent3(k).y(i,Hist_index:end),'FaceColor',c3);
        line5 = histogram(PostfitAgent4(k).y(i,Hist_index:end),'FaceColor',c4);
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
    StateError = [StateError(1,:);1E3*StateError(4,:);StateError(2,:);1E3*StateError(5,:);StateError(3,:);1E3*StateError(6,:)];
    Sigma(k).s = [Sigma(k).s(1,:);1E3*Sigma(k).s(4,:);Sigma(k).s(2,:);1E3*Sigma(k).s(5,:);Sigma(k).s(3,:);1E3*Sigma(k).s(6,:)];
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
    hL = legend('State Error$~~~~~~~~~~~~~~~$','Covariance Envelope');
    % Programatically move the Legend
    newPosition = [.27 .94 0.5 0.05];
    newUnits = 'normalized';
    set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',2);
    % sgt = sgtitle(['State Residual (Agent' num2str(k) ')']);
    % sgt.FontSize = 35;
    
end

%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Plotting the StateError Residuals
close all
Label={'$x$[km]', '$\dot{x}$[m/s]', '$y$[km]', '$\dot{y}$[m/s]', '$z$[km]', '$\dot{z}$[m/s]'};
ColorMatrix = [c1;c6;c3;c4;c5;c2];
figure('Renderer', 'painters', 'Position', [10 10 1000 700])
for k = 1:Num_deputies
    plt   = zeros(6);
    subplot(n/2,2,k)
    plt2 = yline(3,'--','LineWidth', 2, 'Color', 'r');
    hold on
    for i = 1:n
        plt(:,i) = plot(indexY,ConstAnalysis(k).Data(i,1:m),'Color',ColorMatrix(i,:));
    end
    xlabel('Period')
    hold off
    grid on
    clear StateError Sigma2
end
hL = legend([plt(end,1),plt(end,2),plt(end,3),plt2,plt(end,4),plt(end,5),plt(end,6)],...
    {'$x$[km]', '$y$[km]', '$z$[km]','3$\sigma$ Covariance', '$\dot{x}$[m/s]', '$\dot{y}$[m/s]',...
    '$\dot{z}$[m/s]'});
newPosition = [0.68 0.17 0.1 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',2);

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

%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Plotting the trace of the covariance matrices
figure
Label = {'Agent1','Agent2','Agent3','Agent4','Agent5'};
plt   = zeros(5);
for k = 1:Num_deputies
    plt(k) = plot(indexY,TraceVel(k).t(:));
    hold on
    grid on
    legend(plt(1:k),Label(1:k))
end
xlabel('Period')
ylabel('Trace of P')
set(gca, 'YScale', 'log')


