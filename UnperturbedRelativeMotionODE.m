function Xaug_dot = UnperturbedRelativeMotionODE(t,Xaug,mu_Earth,n)
format long 
[l,~] = size(Xaug);
Xaug = reshape(Xaug,n,l/n);

%% Integrating the the target trajectory
Chief           = Xaug(:,1);
r_Chief         = Chief(1:3,:);
v_Chief         = Chief(4:6,:);
rnorm_Chief     = norm(r_Chief);
Xaug_dot(1:3,:) = v_Chief;
Xaug_dot(4:6,:) = -mu_Earth/rnorm_Chief^3*r_Chief;


%% Indexing relevant quatities for the integration
Deputy        = Xaug(:,2:end);
rho           = Deputy(1:3,:);
rho_prime     = Deputy(4:6,:);
r_Deputy      = [rnorm_Chief 0 0]' + rho;
rnorm_Deputy  = vecnorm(r_Deputy); % RIC position of the chasser
TN            = DCM(r_Chief,v_Chief); % DCM's between target's RIC frame and ECI frame

%% Computing the angular velocity and angular acceleration of the target frame
h_vec    = cross(r_Chief,v_Chief);
omega    = TN*( h_vec/rnorm_Chief^2 );
omegadot = TN*( -2*dot(r_Chief,v_Chief)*h_vec/rnorm_Chief^4 );
omega_tilde = [0 -omega(3) omega(2);
    omega(3) 0 -omega(1);
    -omega(2) omega(1) 0];
omegadot_tilde = [0 -omegadot(3) omegadot(2);
    omegadot(3) 0 -omegadot(1);
    -omegadot(2) omegadot(1) 0];

%% Computing the relative acceleration
delf =  -mu_Earth*r_Deputy./(rnorm_Deputy.^3) + mu_Earth/rnorm_Chief^2*[1 0 0]'; % (2-body) gravitatinal

%% Integrating the relative trajectory
Rho_dot(1:3,:) = rho_prime;
Rho_dot(4:6,:) = delf - 2*omega_tilde*rho_prime ...
    - omegadot_tilde*rho - omega_tilde*omega_tilde*rho;

%% Collecting the rates of change
Xaug_dot = [Xaug_dot, Rho_dot];
Xaug_dot = Xaug_dot(:);


end

