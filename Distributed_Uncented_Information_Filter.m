clear all
close all
clc
start_up
format long e
mu_Earth = 3.986004415E5;

%% A) Defining the inital conditions for the spaceraft in the system
for i = 1:5
    
    % % Target initial conditions(Orbital Elements)
    % %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    % % Target initial conditions(Orbital Elements)
    % a1      = 7500;   % 1.5E+4; % Semi-major axis in Km
    % e1      = 0;    % Eccentricity
    % inc1    = 10;     % Inclination in deg0
    % BigOmg1 = 45;     % RAAN in deg
    % LitOmg1 = 10;     % AOP in deg
    % f1      = 0;    % True anomaly in deg
    %
    %
    %
    %
    % % Relative orbital elements of the deputy
    % %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    % dela = 0;                    % Relative semi-major axis in Km
    % deli = i/(4*a1);             % Relative inclination in deg0
    % dele = -1/(2*a1);           % Relative eccentricity, in rad
    % delLitOmg = 0;        % Relative RAAN in rad
    % delBigOmg = 0;               % Relative AOP in rad
    % delf = pi*i;                    % Relative true anomaly in rad
    %
    %
    % % Chaser initial conditions(Orbital Elements)
    % %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    % a2      = a1;          % Semi-major axis in Km
    % e2      = e1 + dele;   % Eccentricity
    % inc2    = inc1;        % Inclination in deg
    % BigOmg2 = BigOmg1;     % RAAN in deg
    % LitOmg2 = LitOmg1;     % AOP in deg
    % f2      = f1;          % True anomaly in deg
    %
    %
    % % Converting angles in radians
    % inc1 = deg2rad(inc1); BigOmg1 = deg2rad(BigOmg1); LitOmg1 = deg2rad(LitOmg1); f1 = deg2rad(f1);
    % COE1 = [a1,e1,inc1,BigOmg1,LitOmg1,f1];
    % [Position_target,Velocity_target]  = COEstoRV(COE1,mu_Earth);
    %
    % inc2 = deg2rad(inc2) + deli; BigOmg2 = deg2rad(BigOmg2); LitOmg2 = deg2rad(LitOmg2)+delLitOmg; f2 = deg2rad(f2);
    % COE2 = [a2,e2,inc2,BigOmg2,LitOmg2,f2];
    % [Position_chaser,Velocity_chaser]  = COEstoRV(COE2,mu_Earth);
    %
    %
    % TN = DCM(Position_target,Velocity_target); rt_norm = norm(Position_target);
    % h_vec = cross(Position_target,Velocity_target); h_norm = norm(h_vec);
    % eh = h_vec/h_norm;
    % N_nudot = h_vec/rt_norm^2;
    % NR_rel = Position_chaser - Position_target; NV_rel2 = Velocity_chaser - Velocity_target;
    % Agents(:,i) = [ TN*NR_rel; TN*(NV_rel2 - cross(N_nudot,NR_rel))];
    
end
a1 = 7500;
Connectivity = [[1, 1, 0, 1, 0]
    [1, 1, 1, 1, 0]
    [0, 1, 1, 1, 1]
    [1, 1, 1, 1, 1]
    [0, 0, 1, 1, 1]];

Agents = [[6.3, -2.3, -2.2, 1.7005e-04   1.4124e-04   -7.3132e-05]
    [-2.1, 4.3, 3.2, 2.7005e-04   1.4124e-06   6.3132e-05]
    [-2.52, -2.67, 8.1, -2.7005e-06   1.4124e-05   3.3132e-07]
    [2.4, -4.6, 4.7, -3.7005e-06   3.4124e-06   7.3132e-06]
    [3.93, 3.7, 2.3, 3.7005e-07   1.6124e-06   8.3132e-07]]';

%% B) Generating the noisy measurements
dt = 5; am = 6;
% Setting the integrator parameters
Period = 2*pi*sqrt(a1^3/mu_Earth);
IntegrationTime = 1*Period;
tvec  = 0:dt:IntegrationTime;
options = odeset('RelTol',2.22045e-13,'AbsTol',2.22045e-30);
[n, Num_agents] = size(Agents);

[~, Xt] = ode113(@(t,X_aug)ClohessyWiltshire_ODE(t, X_aug,a1,mu_Earth,am),tvec,Agents(:),options);
a = 1; b = 6;
for kk = 1:Num_agents
    X_truth(kk).x = Xt(:,a:b)';
    a = a+6;
    b = b+6;
    
    %     plot3(X_truth(kk).x(1,:),X_truth(kk).x(2,:),X_truth(kk).x(3,:));
    %     grid on
    %     hold on
end



Y = nan(4, Num_agents, length(tvec));
sigma1_measurement = 1E-2;
sigma2_measurement = 1E-5;
for i = 1:length(tvec)
    index  = randperm(Num_agents,randi(Num_agents-1));
    target = randi([1 Num_agents]);
    
    if numel(find(index==target)) ~= 0
        while length(index) ~= length(unique(index)) || numel(find(index==target)) ~= 0
            index(index==target) = randi(Num_agents);
            [~, w] = unique( index, 'stable' );
            ind = setdiff( 1:numel(index), w );
            if ~isempty(ind)
                index(ind) = randi(Num_agents);
            end
            
        end
    end
    Xcheck_k = reshape(Xt(i,:),[n, Num_agents]);
    for k = 1:length(index)
        coef = randn(2,1);
        v = [sigma1_measurement; sigma2_measurement].*coef;
        Xi = X_truth(target).x(:,i);
        Xj = X_truth(index(k)).x(:,i);
        Y(:, index(k), i) = [tvec(i); MeasurementFunc(Xi,Xj) + v; target];
    end
end



%% C) Read the next observation: ti,Yi,Ri
m = length(tvec);
Q = zeros(n,n); % diag([(1E-5)^2,(1E-5)^2,(1E-6)^2, (1E-6)^2]); %

%% (1) Calculating all the weights
alpha = 1; % Alpha varies for [1E-4 1]
beta = 2;
L = n; % Number of state to be estimated
kappa = 3-L;
lamda = alpha^2*(L+kappa)-L;
gamma = sqrt(L+lamda);

Wo_m = lamda/(L+lamda);
Wo_c = lamda/(L+lamda)+(1-alpha^2+beta);
Wi = zeros(1,12);
for i = 1:2*L
    Wi(i) = 1/(2*(L+lamda));
end
Wm = [Wo_m Wi];
Wc = [Wo_c Wi];


%% (2) Initializing the UKF

Cov = diag([1 1 1 (1E-3)^2 (1E-3)^2 (1E-3)^2]);
for k = 1:Num_agents
    for kk = 1:Num_agents
        P_post(k,kk).P = Cov; % A priory covariance, at t=0
        X_post(k,kk).x = Agents(:,kk) + Cov*randn(n,1); % Cmplete estimated state, not state deviation, at t=0
    end
end
R = diag([(sigma1_measurement)^2 (sigma2_measurement)^2]); % Measurement covariance

h = waitbar(0,'1','Name','Remaining time...',...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

k = 1;
for i = 1:m % Measurement start at t=0
    
    ti = tvec(i);
    if getappdata(h,'canceling')
        delete(h) % DELETE the waitbar; don't try to CLOSE it.
        wh=findall(0,'tag','TMWWaitbar');
        delete(wh)
        break
    end
    waitbar(i/m,h,sprintf('t = %0.5g',ti))
    
    for k = 1:Num_agents % Referening to the agent doing the tracking
        for kk = 1:Num_agents % Referening the agent being tracked (The kk^th filter in agent k)
            %% (3) Computing Sigma Points (matrix) for t-1
            S_post = sqrtm(P_post(k,kk).P);
            Chi_post = [X_post(k,kk).x, X_post(k,kk).x + gamma*S_post, X_post(k,kk).x - gamma*S_post];
            Num_worker = size(Chi_post,2);
            
            
            %% (4) Propagate Sigma Points through non-linear system
            if ti~=0
                X_aug0 = Chi_post; [am,mm] = size(X_aug0);
                [~, Xnom] = ode113(@(t,X_aug)ClohessyWiltshire_ODE(t, X_aug,a1,mu_Earth,am),[0 dt],X_aug0,options);
                Xnom = reshape(Xnom(end,:),[am,mm]);
            else
                Xnom = Chi_post;
            end
            
            
            %% (5) Perform a state and covariance time update
            X_ap(k,kk).x = zeros(n,1);
            P_ap(k,kk).P = zeros(n,n);
            for j = 1:Num_worker
                Xnomi = Xnom(:,j);
                X_ap(k,kk).x = X_ap(k,kk).x + Wm(j)*Xnomi; % mean of weighted sigma points
            end
            
            for j = 1:Num_worker
                Xnomi = Xnom(:,j);
                P_ap(k,kk).P = P_ap(k,kk).P + Wc(j)*(Xnomi-X_ap(k,kk).x)*(Xnomi-X_ap(k,kk).x)'; % covariance of weighted sigma points
            end
            P_ap(k,kk).P = Q + P_ap(k,kk).P;
            
            %% (6) Re-compute Sigma Points to incorporate effects of process noise
            S_ap = sqrtm(P_ap(k,kk).P);
            Chi_ap(k,kk).x = [X_ap(k,kk).x X_ap(k,kk).x + gamma*S_ap X_ap(k,kk).x - gamma*S_ap];
            
        end
    end
    clear P_post X_post
    for k = 1:Num_agents % Referening to the agent doing the tracking
        for kk = 1:Num_agents % Referening the agent being tracked (The kk^th filter in agent k)
            %% (7) Compute the measurements associated with each Sigma Point vector and their weighted average
            yt_sigma = zeros(2,Num_worker);
            ys_sigma = zeros(2,Num_worker);
            if Y(4, k, i) == kk % Check if agent k took a measurement about agent kk
                for j = 1:Num_worker
                    Xi = Chi_ap(k,kk).x(:,j);
                    Xj = X_ap(k,k).x; % X_truth(k).x(:,i); % 
                    yt_sigma(:,j) = MeasurementFunc(Xi,Xj);
                    
                    Xii = X_ap(k,kk).x;
                    Xjj = Chi_ap(k,k).x(:,j);
                    ys_sigma(:,j) = MeasurementFunc(Xii,Xjj);
                end
            end
            
            
            yt_ap = zeros(2,1);
            ys_ap = zeros(2,1);
            for j = 1:Num_worker
                yi = yt_sigma(:,j);
                yt_ap = yt_ap+Wm(j)*yi;
                
                yii = ys_sigma(:,j);
                ys_ap = ys_ap+Wm(j)*yii;
            end
            
            
            %% (8) Compute the innovation and cross-correlation covariances
            nn = size(yt_ap,1);
            Pt_xy = zeros(n,nn);
            Ps_xy = zeros(n,nn);
            
            for j = 1:Num_worker
                Xnomi = Chi_ap(k,kk).x(:,j);
                yi = yt_sigma(:,j);
                Pt_xy = Pt_xy + Wc(j)*(Xnomi-X_ap(k,kk).x)*(yi-yt_ap)';
                
                
                
                Xnomii = Chi_ap(k,k).x(:,j);
                yii = ys_sigma(:,j);
                Ps_xy = Ps_xy + Wc(j)*(Xnomii-X_ap(k,k).x)*(yii-ys_ap)';
            end
            
            
            %% (9) Compute information vector and the information matrix
            if Y(4, k, i) == kk
                y = Y(2:3, k, i);
            else
                y = zeros(2,1);
            end
            
            % Estimate of the tracking sensor (agent k)
            Hs = 2*(P_ap(k,k).P\Ps_xy)';
            % Computing the information vector and matrix
            H = (P_ap(k,kk).P\Pt_xy)';
            v = y - yt_ap;
            if k == 1
                Info(kk).I_post = H'/(R + Hs*P_ap(k,k).P*Hs')*H;
                info(kk).i_post = H'/(R + Hs*P_ap(k,k).P*Hs')*(v + H*X_ap(k,kk).x);
            else
                Info(kk).I_post = [Info(kk).I_post; H'/(R + Hs*P_ap(k,k).P*Hs')*H];
                info(kk).i_post = [info(kk).i_post; H'/(R + Hs*P_ap(k,k).P*Hs')*(v + H*X_ap(k,kk).x)];
            end
            
        end
    end
    %% (10) Sharing the information gather by each agent via a consensus average
    
    A = Consensus(Connectivity);
    [V,D,W] = eig(A); [~,b] = max(diag(D));
    Weight_Matrix = W(:,b)*V(:,b)';
    Weight_Matrix = kron(Weight_Matrix,eye(n));
    
    for kk = 1:Num_agents
        Info(kk).I_post = Weight_Matrix * Info(kk).I_post;
        info(kk).i_post = Weight_Matrix * info(kk).i_post;
        Info(kk).I_post = permute(reshape(Info(kk).I_post',[n,n,Num_agents]),[2,1,3]);
        info(kk).i_post = reshape(info(kk).i_post,[n,Num_agents]);
    end
    
    %% (11) Performing the measurement update for each agents
    for kk = 1:Num_agents
        PostfitAgent1(kk).y(:,i) = nan(2,1);
        PostfitAgent2(kk).y(:,i) = nan(2,1);
        PostfitAgent3(kk).y(:,i) = nan(2,1);
        PostfitAgent4(kk).y(:,i) = nan(2,1);
        PostfitAgent5(kk).y(:,i) = nan(2,1);
        for k = 1:Num_agents
            omega = ( det(P_ap(k,kk).P\eye(L) + Info(kk).I_post(:,:,k))...
               - det(Info(kk).I_post(:,:,k)) + det(P_ap(k,kk).P\eye(L)))...
                /( 2 * det(P_ap(k,kk).P\eye(L) + Info(kk).I_post(:,:,k)) );
            % Computing the updated state and convariance
            Y_post = omega * P_ap(k,kk).P\eye(L) + (1-omega)*Num_agents * Info(kk).I_post(:,:,k); % information matrix
            P_post(k,kk).P = Y_post\eye(L);
            z_post = omega * P_ap(k,kk).P\X_ap(k,kk).x + (1-omega)*Num_agents * info(kk).i_post(:,k); % information vector
            X_post(k,kk).x = Y_post\z_post;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
            % Computing the residuals
            X_updated(kk).x(:,i) = X_post(k,kk).x;
            if Y(4, k, i)==kk % Check if agent k took a measurement about agent kk
                Xi = X_post(k,kk).x;
                Xj = X_truth(k).x(:,i);
                G = MeasurementFunc(Xi,Xj);
                y = Y(2:3, k, i);
                
                if k==1
                    PostfitAgent1(kk).y(:,i) = y - G;
                elseif k==2
                    PostfitAgent2(kk).y(:,i) = y - G;
                elseif k==3
                    PostfitAgent3(kk).y(:,i) = y - G;
                elseif k==4
                    PostfitAgent4(kk).y(:,i) = y - G;
                elseif k==5
                    PostfitAgent5(kk).y(:,i) = y - G;
                end
                
            end
            % Sampeling the covariance inthe error
            Sigma(kk).s(:,i) = 3*sqrt(diag(P_post(k,kk).P));
        end
    end
    clear P_ap X_ap Chi_ap Xnom_computed
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
close all
c1 = rgb('RosyBrown'); c2 = rgb('Black'); c3 = rgb('Lime');
c4 = rgb('Tomato'); c5 = rgb('DarkBlue');
indexY=1:1:m;
sigma_rho=3*sqrt(R(1,1))*ones(1,m); % Sampeling the covariance inthe error
sigma_rhodot=3*sqrt(R(2,2))*ones(1,m); % Sampeling the covariance inthe error
Label2={'Range' 'Range-Rate'};
Label3={'Range Error','Range-Rate Error'};
Sigma_rho=[sigma_rho;sigma_rhodot];

for kk = 1:Num_agents
    figure('Renderer', 'painters', 'Position', [10 10 900 600])
    for i=1:2
        
        if i==2
            iii=3;
        else
            iii=i;
        end
        subplot(2,2,iii)
        plt1 = plot(indexY,PostfitAgent1(kk).y(i,:),'o','Color', c1);
        hold on
        plt2 = plot(indexY,PostfitAgent2(kk).y(i,:),'o','Color', c2);
        plt3 = plot(indexY,PostfitAgent3(kk).y(i,:),'o','Color', c3);
        plt4 = plot(indexY,PostfitAgent4(kk).y(i,:),'o','Color', c4);
        plt5 = plot(indexY,PostfitAgent5(kk).y(i,:),'o','Color', c5);
        plot(indexY,Sigma_rho(i,:),'r--',indexY,-Sigma_rho(i,:),'r--')
        hold off
        ylabel(Label2(i))
        legend([plt1,plt2,plt3,plt4,plt5],{'Agent1','Agent2','Agent3','Agent4','Agent5'})
        
        
        subplot(2,2,iii+1)
        plt1 = histogram(PostfitAgent1(kk).y(i,:),'FaceColor',c1);
        hold on
        plt2 = histogram(PostfitAgent2(kk).y(i,:),'FaceColor',c2);
        plt3 = histogram(PostfitAgent3(kk).y(i,:),'FaceColor',c3);
        plt4 = histogram(PostfitAgent4(kk).y(i,:),'FaceColor',c4);
        plt5 = histogram(PostfitAgent5(kk).y(i,:),'FaceColor',c5);
        xlabel(Label3(i))
        legend([plt1,plt2,plt3,plt4,plt5],{'Agent1','Agent2','Agent3','Agent4','Agent5'})
        set(gca,'view',[90 -90])
    end
    sgt = sgtitle(['SR-UIF: Measurements Residual With $3\sigma$ Bound (Agent' num2str(kk) ')']);
    sgt.FontSize = 35;
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Plotting the StateError Residuals
Label={'x(t)' 'y(t)' 'z(t)' '$\dot{x}$(t)' '$\dot{y}$(t)' '$\dot{z}$(t)'};
for kk = 1:Num_agents
    StateError = X_updated(kk).x - X_truth(kk).x; % Deviation;
    figure('Renderer', 'painters', 'Position', [10 10 900 600])
    for i=1:n
        subplot(n/2,2,i)
        plot(tvec(1:m),StateError(i,1:m),'b')
        hold on
        plot(tvec(1:m),Sigma(kk).s(i,:),'r--',tvec(1:m),-Sigma(kk).s(i,:),'r--')
        hold off
        ylabel(Label(i))
    end
    sgt = sgtitle(['SR-UIF: State Residual With $3\sigma$ Bound (Agent' num2str(kk) ')']);
    sgt.FontSize = 35;
    
end