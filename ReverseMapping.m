function [InvA] = ReverseMapping(OE,Mu_body)

a = OE(1); theta = OE(2); i = OE(3); q1 = OE(4); q2 = OE(5);

% Defining canstants
p = a*(1 - q1^2 - q2^2);
h = sqrt(Mu_body*p);
R = p/(1 + q1*cos(theta) + q2*sin(theta));
Vr = h/p * (q1*sin(theta) - q2*cos(theta));
Vt = h/p * (1 + q1*cos(theta) + q2*sin(theta));

alpha = a/R;
mu = Vr/Vt;
rho = R/p;
k1 = alpha*(1/rho-1);
k2 = alpha*mu^2/rho;

% Create the reverse mapping [InvA]
InvA = zeros(6,6);

% Nonzero matrix elements of [InvA]
InvA(1,1) = 2*alpha * (2 + 3*k1 + 2*k2);
InvA(1,2) = -2*alpha*mu * ( 1 + 2*k1 + k2);
InvA(1,4) = 2*alpha^2*mu*p / Vt;
InvA(1,5) = 2*a/Vt * ( 1 + 2*k1 + k2 );
InvA(2,2) = 1/R;
InvA(2,3) = cot(i)/R * (cos(theta) + mu*sin(theta));
InvA(2,6) = -sin(theta)*cot(i) / Vt;
InvA(3,3) = (sin(theta) - mu*cos(theta)) / R;
InvA(3,6) = cos(theta) / Vt;
InvA(4,1) = 1/(rho*R) * (3*cos(theta) + 2*mu*sin(theta));
InvA(4,2) = -1/R*( mu^2*sin(theta)/rho + q1*sin(2*theta) - q2*cos(2*theta));
InvA(4,3) = -q2*cot(i)/R * (cos(theta) + mu*sin(theta));
InvA(4,4) = sin(theta)/(rho*Vt);
InvA(4,5) = 1/(rho*Vt) * (2*cos(theta) + mu*sin(theta));
InvA(4,6) = q2*cot(i)*sin(theta)/Vt;
InvA(5,1) = 1/(rho*R) * (3*sin(theta) - 2*mu*cos(theta));
InvA(5,2) = 1/R * (mu^2*cos(theta)/rho + q2*sin(2*theta) + q1*cos(2*theta));
InvA(5,3) = q1*cot(i)/R * (cos(theta) + mu*sin(theta));
InvA(5,4) = -cos(theta)/(rho*Vt);
InvA(5,5) = 1/(rho*Vt) * (2*sin(theta) - mu*cos(theta));
InvA(5,6) = -q1*cot(i)*sin(theta)/Vt;
InvA(6,3) = -(cos(theta) + mu*sin(theta))/(R*sin(i));
InvA(6,6) = sin(theta)/(Vt*sin(i));


end