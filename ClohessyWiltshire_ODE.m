function [Xdot] = ClohessyWiltshire_ODE(t,x,a,mu_Earth,m)
format long
n = sqrt(mu_Earth / a^3);
% State space model
A = [zeros(3) eye(3)
    3*n^2 0 0 0 2*n 0
    0 0 0 -2*n 0 0
    0 0 -n^2 0 0 0];
[l,~] = size(x);
x = reshape(x,m,l/m);
Xdot = A*x;
Xdot = Xdot(:);

end