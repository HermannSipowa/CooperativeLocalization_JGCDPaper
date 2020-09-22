clear all
close all
clc
start_up
format long e

%% A) Initialization for iteration
Connectivity = [[1, 1, 0, 1, 0]
    [1, 1, 1, 1, 0]
    [0, 1, 1, 1, 1]
    [1, 1, 1, 1, 1]
    [0, 0, 1, 1, 1]];
Agents = [[6.3, -2.3, -2.2, 1.7005e-04, 1.4124e-04, -7.3132e-05]
    [-2.1, 4.3, 3.2, 2.7005e-04, 1.4124e-06, 6.3132e-05]
    [-2.52, -2.67, 8.1, -2.7005e-06, 1.4124e-05, 3.3132e-07]
    [2.4, -4.6, 4.7, -3.7005e-06, 3.4124e-06, 7.3132e-06]
    [3.93, 3.7, 2.3, 3.7005e-07, 1.6124e-06, 8.3132e-07]];
Target0 = [-7.0557, -7.3841, 9.9093, 2.7005e-04, 1.4124e-03, 8.3132e-03]';
sigma1_measurement = 1E-2;
sigma2_measurement = 1E-5;


%% B) Generating the noisy measurements
dt = 30; [am,~] = size(Target0);
% Setting the integrator parameters
a1     = 7500;
mu_Earth = 3.986004415E5;
Period = 2*pi*sqrt(a1^3/mu_Earth);
IntegrationTime = 2*Period;
tvec  = 0:dt:IntegrationTime;
options = odeset('RelTol',2.22045e-13,'AbsTol',2.22045e-30);
[~, Xcheck] = ode113(@(t,X_aug)ClohessyWiltshire_ODE(t, X_aug,a1,mu_Earth,am),tvec,Target0,options);


[Num_agents, ~] = size(Agents);
Y = nan(3, Num_agents, length(tvec));
for i = 1:length(tvec)
    index = randi(Num_agents,randi(Num_agents),1);
    for k = 1:length(index)
        coef = randn(2,1);
        v = [sigma1_measurement; sigma2_measurement].*coef;
        Xi = Xcheck(i,:)';
        Xj = Agents(index(k),:)';
        Y(:, index(k), i) = [tvec(i); MeasurementFunc(Xi,Xj) + v];
    end
end



%% C) Read the next observation: ti,Yi,Ri
n = length(Target0); % Size of the state vector
m = length(tvec);
PostfitAgent1 = nan(2,m);
PostfitAgent2 = nan(2,m);
PostfitAgent3 = nan(2,m);
PostfitAgent4 = nan(2,m);
PostfitAgent5 = nan(2,m);
P_post = nan(n, n, Num_agents);
X_post = nan(n, Num_agents);
sigma = nan(n,m);
X_updated = zeros(n,m);
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

for k = 1:Num_agents
    P_post(:,:,k) = diag([1 1 1 (1E-3)^2 (1E-3)^2 (1E-3)^2]); % A priory covariance, at t=0
    error = P_post(:,:,k)*randn(n,1);
    X_post(:,k) = Target0 + error; % Cmplete estimated state, not state deviation, at t=0
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
    Matrix_Pap = nan(n, n, Num_agents);
    Matrix_Xap = nan(n, Num_agents);
    
    for k = 1:Num_agents
        
        %% (3) Computing Sigma Points (matrix) for t-1
        S_post = sqrtm(P_post(:,:,k));
        Chi_post = [X_post(:,k), X_post(:,k)+gamma*S_post, X_post(:,k)-gamma*S_post];
        Num_worker = size(Chi_post,2);
        
        
        %% (4) Propagate Sigma Points through non-linear system
        if ti~=0
            X_aug0 = Chi_post; [am,mm] = size(X_aug0);
            [~, Xnom_computed] = ode113(@(t,X_aug)ClohessyWiltshire_ODE(t, X_aug,a1,mu_Earth,am),[0 dt],X_aug0,options);
            Xnom_computed = reshape(Xnom_computed(end,:),[am,mm]);
        else
            Xnom_computed = Chi_post;
        end
        
        %% (5) Perform a state and covariance time update
        X_ap = zeros(n,1);
        P_ap = zeros(n,n);
        for j = 1:Num_worker
            Xnomi = Xnom_computed(:,j);
            X_ap = X_ap+Wm(j)*Xnomi; % mean of weighted sigma points
        end
        
        for j = 1:Num_worker
            Xnomi = Xnom_computed(:,j);
            P_ap = P_ap+Wc(j)*(Xnomi-X_ap)*(Xnomi-X_ap)'; % covariance of weighted sigma points
        end
        P_ap = Q+P_ap;
        
        %% (6) Re-compute Sigma Points to incorporate effects of process noise
        S_ap = sqrtm(P_ap);
        Chi_ap = [X_ap X_ap+gamma*S_ap X_ap-gamma*S_ap];
        
        
        %% (7) Compute the measurements associated with each Sigma Point vector and their weighted average
        y_computed = zeros(2,Num_worker);
        if isnan(Y(2, k, i)) == 0 % Compute the dynamics if no meassurements initially
            for j = 1:Num_worker
                Xi = Chi_ap(:,j);
                Xj = Agents(k,:)';
                y_computed(:,j) = MeasurementFunc(Xi,Xj);
            end
        end
        
        y_ap = zeros(2,1);
        for j = 1:Num_worker
            yi = y_computed(:,j);
            y_ap = y_ap+Wm(j)*yi;
        end
        
        
        %% (8) Compute the innovation and cross-correlation covariances
        nn = size(y_ap,1);
        P_xy = zeros(n,nn);
        
        for j = 1:Num_worker
            Xnomi = Chi_ap(:,j);
            yi = y_computed(:,j);
            P_xy = P_xy + Wc(j)*(Xnomi-X_ap)*(yi-y_ap)';
        end
        
        
        %% (9) Compute information vector and the information matrix
        if isnan(Y(2, k, i))==0
            y = Y(2:3, k, i);
        else
            y = zeros(2,1);
        end
        
        % Computing the information vector and matrix
        H = (P_ap\P_xy)';
        v = y - y_ap;
        if k == 1
            I_post = H'/R*H;
            i_post = H'/R*(v + H*X_ap);
        else
            I_post = [I_post; H'/R*H];
            i_post = [i_post; H'/R*(v + H*X_ap)];
        end
        Matrix_Pap(:,:,k) = P_ap;
        Matrix_Xap(:,k)   = X_ap;
        
    end
    
    %% (10) Sharing the information gather by each agent via a consensus average
    %     for jj = 1:length(iteration_number)
    %         Weight_Matrix = kron(Consensus(Connectivity),eye(n));
    %         I_post = Weight_Matrix * I_post;
    %         i_post = Weight_Matrix * i_post;
    %     end
    A = Consensus(Connectivity);
    [V,D,W] = eig(A); [~,b] = max(diag(D));
    Weight_Matrix = W(:,b)*V(:,b)';
    Weight_Matrix = kron(Weight_Matrix,eye(n));
    I_post = Weight_Matrix * I_post;
    i_post = Weight_Matrix * i_post;
    I_post = permute(reshape(I_post',[n,n,Num_agents]),[2,1,3]);
    i_post = reshape(i_post,[n,Num_agents]);
    
    
    %% (11) Performing the measurement update for each agents
    for k = 1:Num_agents
        % Computing the updated state and convariance
        Y_post = Matrix_Pap(:,:,k)\eye(L) + Num_agents * I_post(:,:,k); % information matrix
        P_post(:,:,k) = Y_post\eye(L);
        z_post = Matrix_Pap(:,:,k)\Matrix_Xap(:,k) + Num_agents * i_post(:,k); % information vector
        X_post(:,k) = Y_post\z_post;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        % % Computing the residuals
        X_updated(:,i) = X_post(:,k);
        G = nan;
        if isnan(Y(2, k, i)) == 0 % Compute the dynamics if no meassurements initially
            Xi = X_post(:,k);
            Xj = Agents(k,:)';
            G = MeasurementFunc(Xi,Xj);
            y = Y(2:3, k, i);
            if k==1
                PostfitAgent1(:,i) = y - G;
            elseif k==2
                 PostfitAgent2(:,i) = y - G;
            elseif k==3
                 PostfitAgent3(:,i) = y - G;
            elseif k==4
                 PostfitAgent4(:,i) = y - G;
            elseif k==5
                 PostfitAgent5(:,i) = y - G;
            end
            
        end
        % Sampeling the covariance inthe error
        sigma(:,i) = 3*sqrt(diag(P_post(:,:,k)));
    end
end

delete(h) % DELETE the waitbar; don't try to CLOSE it.
wh=findall(0,'tag','TMWWaitbar');
delete(wh)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Plotting the measurements residuals
close all
c1 = rgb('RosyBrown'); c2 = rgb('Black'); c3 = rgb('Lime');
c4 = rgb('Tomato'); c5 = rgb('DarkBlue');
indexY = tvec(1:m)/Period;
sigma_rho=3*sqrt(R(1,1))*ones(1,m); % Sampeling the covariance inthe error
sigma_rhodot=3*sqrt(R(2,2))*ones(1,m); % Sampeling the covariance inthe error
figure('Renderer', 'painters', 'Position', [10 10 900 600])
Label2={'Range' 'Range-Rate'};
Label3={'Range Error','Range-Rate Error'};
Sigma_rho=[sigma_rho;sigma_rhodot];

for i=1:2
    
    if i==2
        iii=3;
    else
        iii=i;
    end
    subplot(2,2,iii)
    plt1 = plot(indexY,PostfitAgent1(i,:),'o','Color', c1);
    hold on
    plt2 = plot(indexY,PostfitAgent2(i,:),'o','Color', c2);
    plt3 = plot(indexY,PostfitAgent3(i,:),'o','Color', c3);
    plt4 = plot(indexY,PostfitAgent4(i,:),'o','Color', c4);
    plt5 = plot(indexY,PostfitAgent5(i,:),'o','Color', c5);
    plot(indexY,Sigma_rho(i,:),'r--',indexY,-Sigma_rho(i,:),'r--')
    hold off
    ylabel(Label2(i))
    xlabel('Period')
    legend([plt1,plt2,plt3,plt4,plt5],{'Agent1','Agent2','Agent3','Agent4','Agent5'})
    
    
    subplot(2,2,iii+1)
    plt1 = histogram(PostfitAgent1(i,:),'FaceColor',c1);
    hold on
    plt2 = histogram(PostfitAgent2(i,:),'FaceColor',c2);
    plt3 = histogram(PostfitAgent3(i,:),'FaceColor',c3);
    plt4 = histogram(PostfitAgent4(i,:),'FaceColor',c4);
    plt5 = histogram(PostfitAgent5(i,:),'FaceColor',c5);
    xlabel(Label3(i))
    legend([plt1,plt2,plt3,plt4,plt5],{'Agent1','Agent2','Agent3','Agent4','Agent5'})
    set(gca,'view',[90 -90])
end
sgt = sgtitle('SR-UIF: Measurements Residual With $3\sigma$ Bound');
sgt.FontSize = 35;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Plotting the StateError Residuals
StateError = X_updated - Xcheck'; % Deviation;
StateError = [StateError(1,:);StateError(4,:);StateError(2,:);StateError(5,:);StateError(3,:);StateError(6,:)];
sigma = [sigma(1,:);sigma(4,:);sigma(2,:);sigma(5,:);sigma(3,:);sigma(6,:)];
Label={'x(t)', '$\dot{x}$(t)', 'y(t)', '$\dot{y}$(t)', 'z(t)', '$\dot{z}$(t)'};
figure('Renderer', 'painters', 'Position', [10 10 900 600])
for i=1:n
    subplot(n/2,2,i)
    plot(indexY,StateError(i,1:m),'b')
    hold on
    plot(indexY,sigma(i,1:m),'r--',indexY,-sigma(i,1:m),'r--')
    hold off
    ylabel(Label(i))
    xlabel('Period')
end
sgt = sgtitle('SR-UIF: State Residual With $3\sigma$ Bound');
sgt.FontSize = 35;
