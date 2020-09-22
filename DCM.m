function [CN] = DCM(r,v)
%--------------------------------------------------------------------------
% [CN] = DCM(Position1,Velocity1)
%
% Relative orientation of the chief relative to the inertial frame
%
% [inputs]: - [r]  : position vector (in Km)
%           - [v]  : Velocity vector (in Km/s)
%
% [outputs]: -[CN] : Direct cosine matrix from ECI to LVLH
%--------------------------------------------------------------------------
x_hat = r/norm(r);
z_hat = cross(r,v)/(norm(cross(r,v)));
y_hat = cross(z_hat,x_hat)/norm(cross(z_hat,x_hat));
CN = [x_hat'; y_hat'; z_hat'];

end