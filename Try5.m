clear all
close all
clc
start_up
format long e
mu_Earth = 3.986004415E5;
c1 = rgb('RosyBrown'); c2 = rgb('Black'); c3 = rgb('Lime');
c4 = rgb('Tomato'); c5 = rgb('DarkBlue');
dt = 10;

%% A) Defining the inital conditions for the spaceraft in the system

% Target initial conditions(Orbital Elements)
%*********************************************%
a1          = 1.5E+4; % Semi-major axis [Km] 7500; % 
e1          = 0;    % Eccentricity
inc1_deg    = 10;   % Inclination [deg]
BigOmg1_deg = 45;   % RAAN [deg]
LitOmg1_deg = 10;   % AOP [deg]
f_deg       = 0;    % True anomaly [deg]

% Converting angles in radians
%******************************%
inc1 = deg2rad(inc1_deg); BigOmg1 = deg2rad(BigOmg1_deg);
LitOmg1 = deg2rad(LitOmg1_deg); f = deg2rad(f_deg);
COE1 = [a1,e1,inc1,BigOmg1,LitOmg1,f];
[Position_target,Velocity_target]  = COEstoRV(COE1,mu_Earth);
X_Chief0 = [Position_target; Velocity_target];
Deputies   = [];
for i = 1:5
    % Relative orbital elements of the deputy
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    dela = 0;                      % Relative semi-major axis [Km]
    dele = -2/(2*a1);              % Relative eccentricity
    deli = 5*i/(4*a1);             % Relative inclination [rad]
    delLitOmg = (-1)^i * 3*i/a1;   % Relative RAAN [rad]
    delBigOmg = (-1)^i * 2*i/a1;   % Relative AOP [rad]
    delf = pi*i/(4*a1);            % Relative true anomaly [rad]
    
    
    % Chaser initial conditions(Orbital Elements)
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    a2          = a1 + dela;       % Semi-major axis [Km]
    e2          = e1 + dele;       % Eccentricity
    inc2_deg    = inc1_deg;        % Inclination  [deg]
    BigOmg2_deg = BigOmg1_deg;     % RAAN [rad]
    LitOmg2_deg = LitOmg1_deg;     % AOP [deg]
    f2_deg      = f_deg;           % True anomaly [deg]
    
    
    inc2 = deg2rad(inc2_deg) + deli; BigOmg2 = deg2rad(BigOmg2_deg)+delBigOmg;
    LitOmg2 = deg2rad(LitOmg2_deg)+delLitOmg; f2 = deg2rad(f2_deg)+delf;
    COE2 = [a2,e2,inc2,BigOmg2,LitOmg2,f2];
    [Position_deputy,Velocity_deputy]  = COEstoRV(COE2,mu_Earth);
    Deputies = [Deputies, [Position_deputy; Velocity_deputy]];
    
end

%% B) Computing all the true trajectories

% Setting the integrator parameters
%***********************************%
[n, Num_deputies] = size(Deputies);
Period = 2*pi*sqrt(a1^3/mu_Earth);
IntegrationTime = 2*Period;
tvec  = 0:dt:IntegrationTime;
options = odeset('RelTol',2.22045e-13,'AbsTol',2.22045e-30);

% Intagrading the chief's and deputies' trajectory
%**************************************************%
[~, X_Chief] = ode113(@(t,X_aug)Unperturbed2Bodyfunc(t, X_aug,mu_Earth,n),tvec,X_Chief0,options);
[~, Xt]      = ode113(@(t,X_aug)Unperturbed2Bodyfunc(t, X_aug,mu_Earth,n),tvec,Deputies,options);
index = 1:6;
for k = 1:Num_deputies
    X_Deputy(k).x = Xt(:,index)';
    index = index + 6;
end

%% C) Visualazing the geomety of the problem

% Trajectory in the ECI frame
%*****************************%
for i = 1
%     plot3(X_Chief(:,1),X_Chief(:,2),X_Chief(:,3))
%     plot3(X_Chief(1,1),X_Chief(1,2),X_Chief(1,3),'o',...
%         'LineWidth',2,...
%         'MarkerEdgeColor','r',...
%         'MarkerFaceColor','r',...
%         'MarkerSize',5);
%     hold on
%     for k = 1:Num_deputies
%     
%         plot3(X_Deputy(k).x(1,:),X_Deputy(k).x(2,:),X_Deputy(k).x(3,:));
%     
%         plot3(X_Deputy(k).x(1,1),X_Deputy(k).x(2,1),X_Deputy(k).x(3,1),'o',...
%             'LineWidth',2,...
%             'MarkerEdgeColor','b',...
%             'MarkerFaceColor','b',...
%             'MarkerSize',5);
%         grid on
%     
%     end
end

% Trajectory in the chief's LVLH frame
%**************************************%
for j =1
%     figure
%     plot3(0,0,0,'bo',...
%         'LineWidth',2,...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor',c2,...
%         'MarkerSize',15);
%     hold on
%     grid on
%     xlabel('X [km]')
%     ylabel('Y [km]')
%     zlabel('Z [km]')
%     for k = 1:Num_deputies
%         for i = 1:length(tvec)
%             Position_target = X_Chief(i,1:3)'; Velocity_target = X_Chief(i,4:6)';
%             Position_deputy = X_Deputy(k).x(1:3,i); Velocity_deputy = X_Deputy(k).x(4:6,i);
%     
%             TN = DCM(Position_target,Velocity_target); rt_norm = norm(Position_target);
%             h_vec = cross(Position_target,Velocity_target); h_norm = norm(h_vec);
%             eh = h_vec/h_norm;
%             N_nudot = h_vec/rt_norm^2;
%             NR_rel = Position_deputy - Position_target; NV_rel2 = Velocity_deputy - Velocity_target;
%             Rel_trajectory(:,i) = [ TN*NR_rel; TN*(NV_rel2 - cross(N_nudot,NR_rel))];
%         end
%     
%         plot3(Rel_trajectory(1,:),Rel_trajectory(2,:),Rel_trajectory(3,:));
%     
%         plot3(Rel_trajectory(1,1),Rel_trajectory(2,1),Rel_trajectory(3,1),'o',...
%             'LineWidth',2,...
%             'MarkerEdgeColor','b',...
%             'MarkerFaceColor','b',...
%             'MarkerSize',5);
%         grid on
%         clear Rel_trajectory
%     end
end

%% D) Generating the noisy measurements

Yrel = nan(4, Num_deputies, length(tvec));
Yabsolute = nan(5, length(tvec));
sigma1_measurement = 1E-2;
sigma2_measurement = 1E-5;

% Absolute measurements from the chief and ground sensors
%*********************************************************%
origin = zeros(6,1);
sensor3 = [1000*ones(3,1); zeros(3,1)];
T_delay = dt; count = 0;
for i = 1:length(tvec)
    if floor(tvec(i)/T_delay)>=count
        index = randi([4 5]);
        sensorID = randi([0 2]);
        coef = randn(2,1);
        v = [sigma1_measurement; sigma2_measurement].*coef;
        Xi = X_Deputy(index).x(:,i);
        if sensorID == 0
            Xj = origin;
        elseif sensorID == 1
            Xj = X_Chief(i,:)';
        elseif sensorID == 2
            Xj = sensor3;
        end
        count = count+1;
        Yabsolute(:, i) = [tvec(i); MeasurementFunc(Xi,Xj) + v; index; sensorID];
    end
end

% Relative measurements between deputies
%****************************************%
for i = 1:length(tvec)
    
    target = randi([1 Num_deputies]);
    index  = randperm(Num_deputies, Num_deputies-1);
    % index  = randperm(Num_agents,randi(Num_agents-1));
    
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
        coef = randn(2,1);
        v = [sigma1_measurement; sigma2_measurement].*coef;
        Xi = X_Deputy(target).x(:,i);
        Xj = X_Deputy(index(k)).x(:,i);
        Yrel(:, index(k), i) = [tvec(i); MeasurementFunc(Xi,Xj) + v; target];
    end
end


%% (1) Initializing the UKF
Cov = diag([1 1 1 (1E-3)^2 (1E-3)^2 (1E-3)^2]);
for k = 1:Num_deputies
    P_post(k).P = Cov;
    X_post(k).x = X_Deputy(k).x(:,1) + Cov*randn(n,1);
end
Q = zeros(n); % diag([(1E-5)^2, (1E-5)^2, (1E-5)^2, (1E-6)^2, (1E-6)^2, (1E-6)^2]);
R = diag([(sigma1_measurement)^2 (sigma2_measurement)^2]); % Measurement covariance


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


h = waitbar(0,'1','Name','Remaining time...',...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

m = length(tvec);
for epoch = 1:m % Measurement start at t=0
    
    ti = tvec(epoch);
    
    if getappdata(h,'canceling')
     delete(h) % DELETE the waitbar; don't try to CLOSE it.
     wh=findall(0,'tag','TMWWaitbar');
     delete(wh)
     break
    end
    waitbar(epoch/m,h,sprintf('t = %0.5g',ti))
    
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
            X_aug0    = Chi_post; [n,mm] = size(X_aug0);
            [~, Xnom] = ode113(@(t,X_aug)Unperturbed2Bodyfunc(t, X_aug,mu_Earth,n),[0 dt],X_aug0,options);
            Xnom = reshape(Xnom(end,:),[n,mm]);
        else
            Xnom = Chi_post;
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
    
    %% (4) Computing estimate using absolute measurements
    %*****************************************************%
    for k = 1:Num_deputies
        Chi_ap = StructChi_ap(k).x; X_ap = StuctX_ap(k).x; P_ap = StructP_ap(k).P;
        Ivector = P_ap\X_ap; Imatrix = inv(P_ap);
        y_computed = zeros(2,Num_worker); y = zeros(2,1);
        if Yabsolute(4, epoch) == k % Checking if there is any absolute sensor measurement
            y = Yabsolute(2:3, epoch);
            if Yabsolute(5, epoch) == 0
                Xj = origin;
            elseif Yabsolute(5, epoch) == 1
                Xj = X_Chief(epoch,:)';
            elseif Yabsolute(5, epoch) == 2
                Xj = sensor3;
            end
            for j = 1:Num_worker
                Xi = Chi_ap(:,j);
                y_computed(:,j) = MeasurementFunc(Xi,Xj);
            end
        end
        
        y_ap = zeros(2,1);
        for j = 1:Num_worker
            yi = y_computed(:,j);
            y_ap = y_ap+Wm(j)*yi;
        end
        
        
        % Compute the innovation and cross-correlation covariances
        %**********************************************************%
        nn = size(y_ap,1);
        P_xy = zeros(n,nn);
        P_yy = zeros(nn,nn);
        
        for j = 1:Num_worker
            Xnomi = Chi_ap(:,j);
            yi = y_computed(:,j);
            P_xy = P_xy + Wc(j)*(Xnomi-X_ap)*(yi-y_ap)';
            P_yy = P_yy + Wc(j)*(yi-y_ap)*(yi-y_ap)';
        end
        
        % Computing the information vector and matrix
        %*********************************************%
        v = y - y_ap;
        
        % P_yy = R+P_yy;
        % K = P_xy/P_yy;
        % StuctX_dist(k).x  = X_ap + K*v;
        % StructP_dist(k).P = P_ap - K*P_yy*K';
        
        
        H = (P_ap\P_xy)';
        Y_post = Imatrix + H'/R*H;            % information matrix
        z_post = Ivector + H'/R*(v + H*X_ap); % information vector
        StuctX_dist(k).x  = Y_post\z_post;
        StructP_dist(k).P = inv(Y_post);
        
        % Re-compute Sigma Points to incorporate effects of process noise
        %*****************************************************************%
        S_dist = sqrtm(StructP_dist(k).P);
        StructChi_dist(k).x = [StuctX_dist(k).x StuctX_dist(k).x ...
                              + gamma*S_dist StuctX_dist(k).x - gamma*S_dist];

        
    end

    %% (4) Message Creation and Sending
    for jj = 1:Num_deputies % Index "jj" is referening to the agent doing the tracking
        Chi_ap_sensor = StructChi_dist(jj).x; 
        X_sensor      = StuctX_dist(jj).x;
        P_ap_sensor   = StructP_dist(jj).P;
        
        for ii = 1:Num_deputies % index "ii" is referening the agent being tracked
            Chi_ap = StructChi_dist(ii).x;
            X_ap   = StuctX_dist(ii).x;
            P_ap   = StructP_dist(ii).P;
            y_target = zeros(2,Num_worker); y = zeros(2,1);
            ys_sigma = zeros(2,Num_worker); ys_ap = zeros(2,1);
            if Yrel(4, jj, epoch) == ii % Checking if agent "jj" took a measurement of agent "ii"
                y = Yrel(2:3, jj, epoch);
                for j = 1:Num_worker
                    Xi = Chi_ap(:,j);
                    Xj = X_sensor;
                    y_target(:,j) = MeasurementFunc(Xi,Xj);
                    
                    Xii = X_ap;
                    Xjj = Chi_ap_sensor(:,j);
                    ys_sigma(:,j) = MeasurementFunc(Xii,Xjj);
                end
            end
            
            y_ap = zeros(2,1);
            for j = 1:Num_worker
                yi = y_target(:,j);
                y_ap = y_ap+Wm(j)*yi;
                
                yii = ys_sigma(:,j);
                ys_ap = ys_ap+Wm(j)*yii;
            end
            
            
            % (4.b) Compute the innovation and cross-correlation covariances
            %***************************************************************%
            nn = size(y_ap,1);
            P_xy = zeros(n,nn);
            Ps_xy = zeros(n,nn);
            
            for j = 1:Num_worker
                Xnomi = Chi_ap(:,j);
                yi = y_target(:,j);
                P_xy = P_xy + Wc(j)*(Xnomi-X_ap)*(yi-y_ap)';
                
                Xnomii = Chi_ap_sensor(:,j);
                yii = ys_sigma(:,j);
                Ps_xy = Ps_xy + Wc(j)*(Xnomii-X_sensor)*(yii-ys_ap)';
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
        CIFused_Info.ivector = StructP_dist(ii).P\StuctX_dist(ii).x;
        CIFused_Info.Imatrix = inv(StructP_dist(ii).P);
        for jj = 1:Num_deputies
            currInfoMessage_ji   = FromNeighbors(jj,ii);
            [CIFused_Info,omega] = CI(CIFused_Info,currInfoMessage_ji);
            omega = omega;
        end
        
        StructP_part(ii).P = inv(CIFused_Info.Imatrix); 
        StructX_part(ii).x = CIFused_Info.Imatrix\CIFused_Info.ivector; 
        
        % (5.2) Re-compute Sigma Points to incorporate effects of process noise
        %***********************************************************************%
        S_part = sqrtm(StructP_part(ii).P);
        StructChi_part(ii).x = [StructX_part(ii).x, StructX_part(ii).x ...
                                + gamma*S_part, StructX_part(ii).x - gamma*S_part];
    end
    
    %% (6) Umpdating the statellite state
    for k = 1:Num_deputies
        Chi_ap = StructChi_part(k).x; X_ap = StructX_part(k).x; P_ap = StructP_part(k).P;
        Ivector = P_ap\X_ap; Imatrix = inv(P_ap);
        y_computed = zeros(2,Num_worker); y = zeros(2,1);
        if Yabsolute(4, epoch) == k % Checking if there is any absolute sensor measurement
            y = Yabsolute(2:3, epoch);
            if Yabsolute(5, epoch) == 0
                Xj = origin;
            elseif Yabsolute(5, epoch) == 1
                Xj = X_Chief(epoch,:)';
            elseif Yabsolute(5, epoch) == 2
                Xj = sensor3;
            end
            for j = 1:Num_worker
                Xi = Chi_ap(:,j);
                y_computed(:,j) = MeasurementFunc(Xi,Xj);
            end
        end
        
        y_ap = zeros(2,1);
        for j = 1:Num_worker
            yi   = y_computed(:,j);
            y_ap = y_ap+Wm(j)*yi;
        end
        
        % Compute the innovation and cross-correlation covariances
        %**********************************************************%
        nn = size(y_ap,1);
        P_xy = zeros(n,nn);
        P_yy=zeros(nn,nn);
        
        for j = 1:Num_worker
            Xnomi = Chi_ap(:,j);
            yi = y_computed(:,j);
            P_xy = P_xy + Wc(j)*(Xnomi-X_ap)*(yi-y_ap)';
            P_yy = P_yy+Wc(j)*(yi-y_ap)*(yi-y_ap)';
        end
        
        % Computing the information vector and matrix
        %*********************************************%
        H = (P_ap\P_xy)';
        v = y - y_ap;
        
        P_yy = R+P_yy;
        K = P_xy/P_yy;
        X_post2 = X_ap + K*v;
        P_post2 = P_ap - K*P_yy*K';
        
        
        Y_post = Imatrix + H'/R*H; % information matrix
        z_post = Ivector + H'/R*(v + H*X_ap); % information vector
        X_post(k).x = Y_post\z_post;
        P_post(k).P = inv(Y_post);
    end
    
    
    %% (7) Computing the measurement residuals
    % Matrix allocation
    %*******************%
    AbsPostfitSensor0.y(:,epoch) = nan(2,1);
    AbsPostfitSensor1.y(:,epoch) = nan(2,1);
    AbsPostfitSensor2.y(:,epoch) = nan(2,1);
    % Computing the residuals for absolute measurements
    %**************************************************%
    for k = 5 % Referening to the agent doing the tracking
        if Yabsolute(4, epoch)==k % Check if agent k took a measurement about agent kk
            if Yabsolute(5,epoch) == 0
                Xj = origin;
            elseif Yabsolute(5,epoch) == 1
                Xj = X_Chief(epoch,:)';
            elseif Yabsolute(5,epoch) == 2
                Xj = sensor3;
            end
            Xi = X_post(k).x;
            G = MeasurementFunc(Xi,Xj);
            y = Yabsolute(2:3,epoch);
            
            if Yabsolute(5,epoch) == 0
                AbsPostfitSensor0.y(:,epoch) = y - G;
            elseif Yabsolute(5,epoch) == 1
                AbsPostfitSensor1.y(:,epoch) = y - G;
            elseif Yabsolute(5,epoch) == 2
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
        PostfitAgent1(k).y(:,epoch) = nan(2,1);
        PostfitAgent2(k).y(:,epoch) = nan(2,1);
        PostfitAgent3(k).y(:,epoch) = nan(2,1);
        PostfitAgent4(k).y(:,epoch) = nan(2,1);
        PostfitAgent5(k).y(:,epoch) = nan(2,1);
        
        % Computing the residuals
        %*************************%
        for j = 1:Num_deputies % Referening to the agent doing the tracking
            if j==k
                continue
            end
            if Yrel(4,j,epoch)==k % Check if agent k took a measurement about agent kk
                Xi = X_post(k).x;
                Xj = X_Deputy(j).x(:,epoch); % X_post(j).x;
                G = MeasurementFunc(Xi,Xj);
                y = Yrel(2:3,j,epoch);
                
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
    end
    
    clear P_ap X_ap Chi_ap CIFused_Info FromNeighbors
end

delete(h) % DELETE the waitbar; don't try to CLOSE it.
wh=findall(0,'tag','TMWWaitbar');
delete(wh)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%
% Plotting the measurements residuals
% close all
indexY = tvec(1:m)/Period;
sigma_rho=3*sqrt(R(1,1))*ones(1,m); % Sampeling the covariance inthe error
sigma_rhodot=3*sqrt(R(2,2))*ones(1,m); % Sampeling the covariance inthe error
Label2={'Range' 'Range-Rate'};
Label3={'Range Error','Range-Rate Error'};
Sigma_rho=[sigma_rho;sigma_rhodot];

for k = 5
    figure('Renderer', 'painters', 'Position', [10 10 900 600])
    for i=1:2
        
        if i==2
            iii=3;
        else
            iii=i;
        end
        subplot(2,2,iii)
        plt1 = plot(indexY,AbsPostfitSensor0.y(i,:),'o','Color', c1);
        hold on
        plt2 = plot(indexY,AbsPostfitSensor1.y(i,:),'o','Color', c2);
        plt3 = plot(indexY,AbsPostfitSensor2.y(i,:),'o','Color', c3);
        plot(indexY,Sigma_rho(i,:),'r--',indexY,-Sigma_rho(i,:),'r--')
        hold off
        ylabel(Label2(i))
        xlabel('Period')
        legend([plt1,plt2,plt3],{'Sensor0','Sensor1','Sensor2'})
        
        
        subplot(2,2,iii+1)
        plt1 = histogram(AbsPostfitSensor0.y(i,:),'FaceColor',c1);
        hold on
        plt2 = histogram(AbsPostfitSensor1.y(i,:),'FaceColor',c2);
        plt3 = histogram(AbsPostfitSensor2.y(i,:),'FaceColor',c3);
        xlabel(Label3(i))
        legend([plt1,plt2,plt3],{'Sensor0','Sensor1','Sensor2'})
        set(gca,'view',[90 -90])
    end
    sgt = sgtitle(['Absolute Measurements Residual (Agent' num2str(k) ')']);
    sgt.FontSize = 35;
    
end

for k = 1:Num_deputies
    figure('Renderer', 'painters', 'Position', [10 10 900 600])
    for i=1:2
        
        if i==2
            iii=3;
        else
            iii=i;
        end
        subplot(2,2,iii)
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
        legend([plt1,plt2,plt3,plt4,plt5],{'Agent1','Agent2','Agent3','Agent4','Agent5'})
        
        
        subplot(2,2,iii+1)
        plt1 = histogram(PostfitAgent1(k).y(i,:),'FaceColor',c1);
        hold on
        plt2 = histogram(PostfitAgent2(k).y(i,:),'FaceColor',c2);
        plt3 = histogram(PostfitAgent3(k).y(i,:),'FaceColor',c3);
        plt4 = histogram(PostfitAgent4(k).y(i,:),'FaceColor',c4);
        plt5 = histogram(PostfitAgent5(k).y(i,:),'FaceColor',c5);
        xlabel(Label3(i))
        legend([plt1,plt2,plt3,plt4,plt5],{'Agent1','Agent2','Agent3','Agent4','Agent5'})
        set(gca,'view',[90 -90])
    end
    sgt = sgtitle(['Relative Measurements Residual (Agent' num2str(k) ')']);
    sgt.FontSize = 35;
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Plotting the StateError Residuals
Label={'x(t)', '$\dot{x}$(t)', 'y(t)', '$\dot{y}$(t)', 'z(t)', '$\dot{z}$(t)'};
for k = 1:Num_deputies
    StateError = X_updated(k).x - X_Deputy(k).x(:,1:m);
    StateError = [StateError(1,:);StateError(4,:);StateError(2,:);StateError(5,:);StateError(3,:);StateError(6,:)];
    Sigma(k).s = [Sigma(k).s(1,:);Sigma(k).s(4,:);Sigma(k).s(2,:);Sigma(k).s(5,:);Sigma(k).s(3,:);Sigma(k).s(6,:)];
    figure('Renderer', 'painters', 'Position', [10 10 900 600])
    for i=1:n
        subplot(n/2,2,i)
        plot(indexY,StateError(i,1:m),'b')
        hold on
        plot(indexY,Sigma(k).s(i,:),'g--',indexY,-Sigma(k).s(i,:),'r--')
        hold off
        ylabel(Label(i))
        xlabel('Period')
    end
    sgt = sgtitle(['State Residual (Agent' num2str(k) ')']);
    sgt.FontSize = 35;
    
end

