function [CIFused,omega,fmin] = CI(Agent_i,Agent_j)

format long
I1 = Agent_i.Imatrix; i1 = Agent_i.ivector;
I2 = Agent_j.Imatrix; i2 = Agent_j.ivector;
f = @(omega) det(inv(omega * I1 + (1-omega) * I2) / det(inv(I1)) );
options = optimset('TolFun',2.22045e-14,'TolX',2.22045e-50,'Display','off');
[omega,fmin] = fminbnd(f,0,1,options);
CIFused.Imatrix = omega * I1 + (1-omega) * I2;
CIFused.ivector = omega * i1 + (1-omega) * i2;


% if norm(i2)==0
%     CIFused.Imatrix = I1;
%     CIFused.ivector = i1;
%     omega = 1;
%     fmin = 0;
% else
%     f = @(omega) det(inv((omega*I1+(1-omega)*I2)))/det(inv(I1));
%     options = optimset('TolFun',2.22045e-14,'TolX',2.22045e-50,'Display','off');
%     [omega,fmin] = fminbnd(f,0,1,options);
%     % [x,fval]     = fminsearch(f2,0,options);
%     % omega = (det(I1+I2)-det(I2)+det(I1))/(2*det(I1+I2));
%     % figure
%     % fplot(f,[0 20],'b','Linewidth',2)
%     % set(gca, 'YScale', 'log')
%     
%     CIFused.Imatrix = omega * I1 + (1-omega) * I2;
%     CIFused.ivector = omega * i1 + (1-omega) * i2;
% end


end