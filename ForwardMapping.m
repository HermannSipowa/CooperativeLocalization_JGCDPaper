function [A] = ForwardMapping(OE,Mu_body)

a = OE(1); theta = OE(2); i = OE(3); q1 = OE(4); q2 = OE(5);

% Defining canstants
p = a*(1 - q1^2 - q2^2);
h = sqrt(Mu_body*p);
R = p/(1 + q1*cos(theta) + q2*sin(theta));
Vr = h/p * (q1*sin(theta) - q2*cos(theta));
Vt = h/p * (1 + q1*cos(theta) + q2*sin(theta));


% Create the forward mapping [A]
A = zeros(6,6);

% Nonzero matrix elements of [A]
A(1,1) = R/a;
A(1,2) = Vr/Vt*R;
A(1,4) = -R/p * ( 2*a*q1 + R*cos(theta) );
A(1,5) = -R/p * ( 2*a*q2 + R*sin(theta) );
A(2,2) = R;
A(2,6) = R*cos(i);
A(3,3) = R*sin(theta);
A(3,6) = -R*cos(theta)*sin(i);
A(4,1) = -Vr/(2*a);
A(4,2) = ( 1/R - 1/p )*h;
A(4,4) = ( Vr*a*q1 + h*sin(theta) )/p;
A(4,5) = ( Vr*a*q2 - h*cos(theta) )/p;
A(5,1) = -3*Vt/(2*a);
A(5,2) = -Vr;
A(5,4) = ( 3*Vt*a*q1 + 2*h*cos(theta) )/p;
A(5,5) = ( 3*Vt*a*q2 + 2*h*sin(theta) )/p;
A(5,6) = Vr*cos(i);
A(6,3) = Vt*cos(theta) + Vr*sin(theta);
A(6,6) = ( Vt*sin(theta) -Vr*cos(theta) )*sin(i);

end