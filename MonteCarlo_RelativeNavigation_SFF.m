clear variables
close all
delete MonteCarloTraceData.mat
delete Agent1_ConstAnalysis.mat
delete Agent2_ConstAnalysis.mat
delete Agent3_ConstAnalysis.mat
delete Agent4_ConstAnalysis.mat
delete Agent5_ConstAnalysis.mat
clc
system('caffeinate -dims &');
format long e
mu_Earth = 3.986004415E5;
c1 = rgb('RosyBrown'); c2 = rgb('Black'); c3 = rgb('Lime');
c4 = rgb('Tomato'); c5 = rgb('DarkBlue'); c6 = rgb('DarkTurquoise');
dt = 60;
MC_runs = 1000; % Number of Monte Carlo simulations

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
indexY            = tvec/Period;
save('indexYMonteCarlo.mat','indexY')
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


%% Monte Carlo
% progressbar('Monte Carlo Trials','Simulation') % Init 2 bars
for Monte = 1:MC_runs
    
    progressbar([],0) % Reset 2nd bar
    
    P_post = struct('P', cell(1, Num_deputies));
    X_post = struct('x', cell(1, Num_deputies));
    StuctX_ap = struct('x', cell(1, Num_deputies));
    StructP_ap = struct('P', cell(1, Num_deputies));
    StructChi_ap = struct('P', cell(1, Num_deputies));
    X_updated = struct('x', cell(1, Num_deputies));
    Sigma     = struct('s', cell(1, Num_deputies));
    TracePos  = struct('t', cell(1, Num_deputies));
    TraceVel  = struct('t', cell(1, Num_deputies));
    %% D) Generating the noisy measurements
    p = 4;
    n = 6;
    Yrel = nan(p+2, Num_deputies, length(tvec));
    Yabsolute = nan(p+3, length(tvec));
    sigma1_measurement = 1E-3;
    sigma2_measurement = 1E-6;
    sigma3_measurement = 1E-3;
    sigma4_measurement = 1E-3;
    
    % Relative measurements between deputies
    %****************************************%
    for i = 1:length(tvec)
        target = randi([1 Num_deputies]);
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
        FromNeighbors = struct('ivector', cell(Num_deputies, Num_deputies), 'Imatrix',cell(Num_deputies, Num_deputies));
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
        
        %% (7) Computing the state residuals
        for k = 1:Num_deputies
            % Storing estimated state
            %*************************%
            X_updated(k).x(:,epoch) = X_post(k).x;
            % Sampeling the covariance inthe error
            %**************************************%
            Sigma(k).s(:,epoch) = 3*sqrt(diag(P_post(k).P));
            TracePos(k).t(epoch) = trace(P_post(k).P(1:3,1:3));
            TraceVel(k).t(epoch) = trace(P_post(k).P(4:6,4:6));
        end
        
        clear P_ap X_ap Chi_ap CIFused_Info FromNeighbors
    end
    
    % Storing data for the consistency monte carlo simulation
    ConstAnalysis = struct('Data', cell(1, Num_deputies));
    for k = 1:Num_deputies
        StateError = X_updated(k).x - X_Deputy(k).x(:,1:m);
        StateError = [StateError(1,:); StateError(2,:); StateError(3,:); StateError(4,:); StateError(5,:); StateError(6,:)];
        ConstAnalysis(k).Data = StateError./Sigma(k).s;
    end
    
    %% Storing the trace of the convariance matrice for every element
    if Monte == 1
        % Storing the trace of the covariance at each MC iteration
        PosData_1 = [TracePos(1).t(:)';TracePos(2).t(:)';TracePos(3).t(:)';...
            TracePos(4).t(:)';TracePos(5).t(:)'];
        save('MonteCarloTraceData.mat','PosData_1')
        TraceMat = matfile('MonteCarloTraceData.mat','Writable',true);
        TraceMat.VelData_1 = [TraceVel(1).t(:)';TraceVel(2).t(:)';TraceVel(3).t(:)';...
            TraceVel(4).t(:)';TraceVel(5).t(:)'];
        
        % Storing the normalized error for each agent at each MC iteration
        Run_1 =  ConstAnalysis(1).Data;
        save('Agent1_ConstAnalysis.mat','Run_1')
        Agent1Const = matfile('Agent1_ConstAnalysis.mat','Writable',true);
        
        Run_1 = ConstAnalysis(2).Data;
        save('Agent2_ConstAnalysis.mat','Run_1')
        Agent2Const = matfile('Agent2_ConstAnalysis.mat','Writable',true);
        
        Run_1 = ConstAnalysis(3).Data;
        save('Agent3_ConstAnalysis.mat','Run_1')
        Agent3Const = matfile('Agent3_ConstAnalysis.mat','Writable',true);
        
        Run_1 = ConstAnalysis(4).Data;
        save('Agent4_ConstAnalysis.mat','Run_1')
        Agent4Const = matfile('Agent4_ConstAnalysis.mat','Writable',true);
        
        Run_1 = ConstAnalysis(5).Data;
        save('Agent5_ConstAnalysis.mat','Run_1')
        Agent5Const = matfile('Agent5_ConstAnalysis.mat','Writable',true);
        
    else
        Posfield = strcat('PosData_',num2str(Monte));
        TraceMat.(Posfield) = [TracePos(1).t(:)';TracePos(2).t(:)';TracePos(3).t(:)';...
            TracePos(4).t(:)';TracePos(5).t(:)'];
        Velfield = strcat('VelData_',num2str(Monte));
        TraceMat.(Velfield) = [TraceVel(1).t(:)';TraceVel(2).t(:)';TraceVel(3).t(:)';...
            TraceVel(4).t(:)';TraceVel(5).t(:)'];
        
        Constfield = strcat('Run_',num2str(Monte));
        Agent1Const.(Constfield) = ConstAnalysis(1).Data;
        Agent2Const.(Constfield) = ConstAnalysis(2).Data;
        Agent3Const.(Constfield) = ConstAnalysis(3).Data;
        Agent4Const.(Constfield) = ConstAnalysis(4).Data;
        Agent5Const.(Constfield) = ConstAnalysis(5).Data;
        
    end
    clear Yrel P_post X_post Cov R Q X_updated StateError Sigma ConstAnalysis TracePos TraceVel
    progressbar(Monte/MC_runs) % Update 1st bar
    
end % Ending the Monte Carlo loop
system('killall caffeinate');

