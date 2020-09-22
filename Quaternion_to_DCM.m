function C = Quaternion_to_DCM(x)

if (size(x,1) ~= 4)
    error('Quaternion parameter set must be 4-by-n matrix.')
end
if (abs(norm(x)-1)> 0.000001)
    error('Quaternions must have unit norm.')
end
C = [x(1)^2+x(2)^2-x(3)^2-x(4)^2 2*(x(2)*x(3)+x(1)*x(4)) 2*(x(2)*x(4)-x(1)*x(3));...
    2*(x(2)*x(3)-x(1)*x(4)) x(1)^2-x(2)^2+x(3)^2-x(4)^2 2*(x(3)*x(4)+x(1)*x(2));...
    2*(x(2)*x(4)+x(1)*x(3)) 2*(x(3)*x(4)-x(1)*x(2)) x(1)^2-x(2)^2-x(3)^2+x(4)^2];

end