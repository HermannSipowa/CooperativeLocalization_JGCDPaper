function Xaug_dot = RelativeMotionODE_cannonball(t, Xaug, Target, Chasser)
format long
global mu_Earth

%% Collecting relevant quantities
Target_States = Xaug(1:6,:); 
rho_aug  = Xaug(7:12,:);

r_target = Target_States(1:3,:); v_target = Target_States(4:6,:);
rt_norm = norm(r_target);
rho = rho_aug(1:3,:); rho_prime = rho_aug(4:6,:);
r_c = [rt_norm+rho(1) rho(2) rho(3)]'; rc_norm = norm(r_c); % RIC position of the chasser
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame
r_chasser = TN\r_c; % ECI position of the chasser


%% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, V_sun, beta_chief] = F_CanonBall(t,r_target,Target); % SRP force on the Target
U_eci_Chasser = F_CanonBall(t,r_chasser,Chasser); % SRP force on the Chasser

U_eci_Chasser = 0*U_eci_Chasser; U_eci_Target = 0*U_eci_Target;
%% Computing the angular velocity and angular acceleration of the target frame
h_vec = cross(r_target,v_target); h_norm = norm(h_vec);  u_acc = U_eci_Target;
% U_eci_Target = u_acc;
eh = h_vec/h_norm; h_dot = cross(r_target,u_acc); r_rel_norm = norm(Xrel);
eh_dot = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
V_rel = V_sun - v_target;
u_acc_dot = 0*beta_chief*(V_rel/r_rel_norm^3 - 3*dot(Xrel,V_rel)*Xrel/r_rel_norm^5);

    
Omega = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
Omegadot = TN*( h_dot/rt_norm^2 ...
                -2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(u_acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );


%% Computing the relative acceleration
delf     =  -mu_Earth/rc_norm^3*r_c + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
DelU_srp = TN*(U_eci_Chasser - U_eci_Target); % Delta-SRP

%% Integrating the the target trajectory
Xdot(1:3,:) = v_target;
Xdot(4:6,:) = -mu_Earth/rt_norm^3*r_target + U_eci_Target;

%% Integrating the relative trajectory
Rho_dot(1:3,:) = rho_prime;
Rho_dot(4:6,:) = DelU_srp + delf - 2*cross(Omega,rho_prime) ...
                - cross(Omegadot,rho) - cross(Omega,cross(Omega,rho));

%% Collecting the rates of change
Xaug_dot = [Xdot; Rho_dot];


end