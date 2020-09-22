function [OBE] = RVtoCOEs(r,v,mu)
%--------------------------------------------------------------------------
% [OBE]=RVtoCOEs(r,v,nu)
%
% This function coverts cartesian corrdinates into classical orbital
% elements (in radian)
% [inputs]: - [r]  : position vector (in Km)
%           - [v]  : Velocity vector (in Km/s)
%           - [nu] : Body's gravitational constant
%
% [outputs]: -[OBE] : Classical Obital Elements (angles in radian)
%--------------------------------------------------------------------------

k = [0 0 1];    % normal vector in the k-direction
h = cross(r,v); % Anglar Momentum Vector
n = cross(k,h); % Node Vector
e_vec = (1/mu)*((norm(v)^2-mu/norm(r))*r-dot(r,v)*v); % Eccentricity Vector


P = norm(h)^2/mu; % Semi-Latus Rectum

e = norm(e_vec); % Eccentricity

i = acos(h(3)/norm(h)); % Inclination

Big_Omega = acos(n(1)/norm(n)); % Longitude of ascending node
if n(2)<0
    Big_Omega = 2*pi-Big_Omega;
end

Little_Omega = acos(dot(n,e_vec)/(norm(n)*norm(e_vec))); % Argument of periapsis
if e_vec(3)<0
    Little_Omega = 2*pi-Little_Omega;
end

Nu=acos(dot(e_vec,r)/(norm(e_vec)*norm(r))); % True Anomaly
if dot(r,v)<0
    Nu = 2*pi-Nu;
end

a = P/(1-norm(e_vec)^2); % Semi-major Axis

OBE = [a; e; i; Big_Omega; Little_Omega; Nu];
end