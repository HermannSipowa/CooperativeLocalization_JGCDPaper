function y = MeasurementFunc(Xi,Xj)

% Extracting the relative vectors 
X_agent_i = Xi(1:3); V_agent_i = Xi(4:6);
X_agent_j = Xj(1:3); V_agent_j = Xj(4:6);
Xrel =  X_agent_j - X_agent_i;
Vrel =  V_agent_j - V_agent_i;

% Computing the relative range and range-rate
rho = norm(Xrel);
rhodot = dot(Xrel,Vrel)/rho;

% Comuting the relative azimuth and elevation
b1 = X_agent_j/norm(X_agent_j);
b3 = cross(X_agent_j,V_agent_j)/norm(cross(X_agent_j,V_agent_j));
b2 = cross(b3,b1);
CR = [b1';b2'; b3'];
xyz = CR*Xrel;
El = asin(xyz(3)/norm(xyz));
Az = atan2(xyz(1),xyz(2));

y = [rho,rhodot, El, Az]';

end