function [Xdot] = Unperturbed2Bodyfunc(t,X,mu_Earth,m)
format long
[l,~] = size(X);
X = reshape(X,m,l/m);
Xdot(1:3,:) = X(4:6,:);
Xdot(4:6,:) = -mu_Earth./vecnorm(X(1:3,:),2,1).^3.*X(1:3,:);
Xdot = Xdot(:);
end