delete(h) % DELETE the waitbar; don't try to CLOSE it.
wh=findall(0,'tag','TMWWaitbar');
delete(wh)


% Absolute measurements from the chief and ground sensors
%*********************************************************%
% origin = zeros(6,1);
% sensor3 = [1000*ones(3,1); zeros(3,1)];
% T_delay = dt; count = 0;
for i = 1:length(tvec)
    %     if floor(tvec(i)/T_delay)>=count
    %         index = randi([4 5]);
    %         sensorID = randi([0 2]);
    %         coef = randn(2,1);
    %         v = [sigma1_measurement; sigma2_measurement].*coef;
    %         Xi = X_Deputy(index).x(:,i);
    %         if sensorID == 0
    %             Xj = origin;
    %         elseif sensorID == 1
    %             Xj = X_Chief(i,:)';
    %         elseif sensorID == 2
    %             Xj = sensor3;
    %         end
    %         count = count+1;
    %         Yabsolute(:, i) = [tvec(i); MeasurementFunc(Xi,Xj) + v; index; sensorID];
    %     end
end

% (4) Computing estimate using absolute measurements
%*****************************************************%
for k = 1:Num_deputies
    % Chi_ap = StructChi_ap(k).x; X_ap = StuctX_ap(k).x; P_ap = StructP_ap(k).P;
    % Ivector = P_ap\X_ap; Imatrix = inv(P_ap);
    % y_computed = zeros(2,Num_worker); y = zeros(2,1);
    % if Yabsolute(4, epoch) == k % Checking if there is any absolute sensor measurement
    %     y = Yabsolute(2:3, epoch);
    %     if Yabsolute(5, epoch) == 0
    %         Xj = origin;
    %     elseif Yabsolute(5, epoch) == 1
    %         Xj = X_Chief(epoch,:)';
    %     elseif Yabsolute(5, epoch) == 2
    %         Xj = sensor3;
    %     end
    %     for j = 1:Num_worker
    %         Xi = Chi_ap(:,j);
    %         y_computed(:,j) = MeasurementFunc(Xi,Xj);
    %     end
    % end
    %
    % y_ap = zeros(2,1);
    % for j = 1:Num_worker
    %     yi = y_computed(:,j);
    %     y_ap = y_ap+Wm(j)*yi;
    % end
    %
    %
    % % Compute the innovation and cross-correlation covariances
    % %**********************************************************%
    % nn = size(y_ap,1);
    % P_xy = zeros(n,nn);
    % P_yy = zeros(nn,nn);
    %
    % for j = 1:Num_worker
    %     Xnomi = Chi_ap(:,j);
    %     yi = y_computed(:,j);
    %     P_xy = P_xy + Wc(j)*(Xnomi-X_ap)*(yi-y_ap)';
    %     P_yy = P_yy + Wc(j)*(yi-y_ap)*(yi-y_ap)';
    % end
    %
    % % Computing the information vector and matrix
    % %*********************************************%
    % v = y - y_ap;
    %
    % % P_yy = R+P_yy;
    % % K = P_xy/P_yy;
    % % StuctX_dist(k).x  = X_ap + K*v;
    % % StructP_dist(k).P = P_ap - K*P_yy*K';
    %
    %
    % H = (P_ap\P_xy)';
    % Y_post = Imatrix + H'/R*H;            % information matrix
    % z_post = Ivector + H'/R*(v + H*X_ap); % information vector
    % StuctX_dist(k).x  = Y_post\z_post;
    % StructP_dist(k).P = inv(Y_post);
    %
    % % Re-compute Sigma Points to incorporate effects of process noise
    % %*****************************************************************%
    % S_dist = sqrtm(StructP_dist(k).P);
    % StructChi_dist(k).x = [StuctX_dist(k).x StuctX_dist(k).x ...
    %     + gamma*S_dist StuctX_dist(k).x - gamma*S_dist];
end



% (5.2) Re-compute Sigma Points to incorporate effects of process noise
%***********************************************************************%
% S_part = sqrtm(StructP_part(ii).P);
% StructChi_part(ii).x = [StructX_part(ii).x, StructX_part(ii).x ...
%     + gamma*S_part, StructX_part(ii).x - gamma*S_part];


% (6) Umpdating the statellite state
for k = 1:Num_deputies
    % Chi_ap = StructChi_part(k).x; X_ap = StructX_part(k).x; P_ap = StructP_part(k).P;
    % Ivector = P_ap\X_ap; Imatrix = inv(P_ap);
    % y_computed = zeros(2,Num_worker); y = zeros(2,1);
    % if Yabsolute(4, epoch) == k % Checking if there is any absolute sensor measurement
    %     y = Yabsolute(2:3, epoch);
    %     if Yabsolute(5, epoch) == 0
    %         Xj = origin;
    %     elseif Yabsolute(5, epoch) == 1
    %         Xj = X_Chief(epoch,:)';
    %     elseif Yabsolute(5, epoch) == 2
    %         Xj = sensor3;
    %     end
    %     for j = 1:Num_worker
    %         Xi = Chi_ap(:,j);
    %         y_computed(:,j) = MeasurementFunc(Xi,Xj);
    %     end
    % end
    %
    % y_ap = zeros(2,1);
    % for j = 1:Num_worker
    %     yi   = y_computed(:,j);
    %     y_ap = y_ap+Wm(j)*yi;
    % end
    %
    % % Compute the innovation and cross-correlation covariances
    % %**********************************************************%
    % nn = size(y_ap,1);
    % P_xy = zeros(n,nn);
    % P_yy=zeros(nn,nn);
    %
    % for j = 1:Num_worker
    %     Xnomi = Chi_ap(:,j);
    %     yi = y_computed(:,j);
    %     P_xy = P_xy + Wc(j)*(Xnomi-X_ap)*(yi-y_ap)';
    %     P_yy = P_yy+Wc(j)*(yi-y_ap)*(yi-y_ap)';
    % end
    %
    % % Computing the information vector and matrix
    % %*********************************************%
    % H = (P_ap\P_xy)';
    % v = y - y_ap;
    %
    % P_yy = R+P_yy;
    % K = P_xy/P_yy;
    % X_post2 = X_ap + K*v;
    % P_post2 = P_ap - K*P_yy*K';
    %
    %
    % Y_post = Imatrix + H'/R*H; % information matrix
    % z_post = Ivector + H'/R*(v + H*X_ap); % information vector
    % X_post(k).x = Y_post\z_post;
    % P_post(k).P = inv(Y_post);
end


% Matrix allocation
%*******************%
% AbsPostfitSensor0.y(:,epoch) = nan(2,1);
% AbsPostfitSensor1.y(:,epoch) = nan(2,1);
% AbsPostfitSensor2.y(:,epoch) = nan(2,1);
% Computing the residuals for absolute measurements
%**************************************************%
for k = 5 % Referening to the agent doing the tracking
    % if Yabsolute(4, epoch)==k % Check if agent k took a measurement about agent kk
    %     if Yabsolute(5,epoch) == 0
    %         Xj = origin;
    %     elseif Yabsolute(5,epoch) == 1
    %         Xj = X_Chief(epoch,:)';
    %     elseif Yabsolute(5,epoch) == 2
    %         Xj = sensor3;
    %     end
    %     Xi = X_post(k).x;
    %     G = MeasurementFunc(Xi,Xj);
    %     y = Yabsolute(2:3,epoch);
    %
    %     if Yabsolute(5,epoch) == 0
    %         AbsPostfitSensor0.y(:,epoch) = y - G;
    %     elseif Yabsolute(5,epoch) == 1
    %         AbsPostfitSensor1.y(:,epoch) = y - G;
    %     elseif Yabsolute(5,epoch) == 2
    %         AbsPostfitSensor2.y(:,epoch) = y - G;
    %     end
    % end
end

for k = 5
    % figure('Renderer', 'painters', 'Position', [10 10 900 600])
    % for i=1:2
    %
    %     if i==2
    %         iii=3;
    %     else
    %         iii=i;
    %     end
    %     subplot(2,2,iii)
    %     plt1 = plot(indexY,AbsPostfitSensor0.y(i,:),'o','Color', c1);
    %     hold on
    %     plt2 = plot(indexY,AbsPostfitSensor1.y(i,:),'o','Color', c2);
    %     plt3 = plot(indexY,AbsPostfitSensor2.y(i,:),'o','Color', c3);
    %     plot(indexY,Sigma_rho(i,:),'r--',indexY,-Sigma_rho(i,:),'r--')
    %     hold off
    %     ylabel(Label2(i))
    %     xlabel('Period')
    %     legend([plt1,plt2,plt3],{'Sensor0','Sensor1','Sensor2'})
    %
    %
    %     subplot(2,2,iii+1)
    %     plt1 = histogram(AbsPostfitSensor0.y(i,:),'FaceColor',c1);
    %     hold on
    %     plt2 = histogram(AbsPostfitSensor1.y(i,:),'FaceColor',c2);
    %     plt3 = histogram(AbsPostfitSensor2.y(i,:),'FaceColor',c3);
    %     xlabel(Label3(i))
    %     legend([plt1,plt2,plt3],{'Sensor0','Sensor1','Sensor2'})
    %     set(gca,'view',[90 -90])
    % end
    % sgt = sgtitle(['Absolute Measurements Residual (Agent' num2str(k) ')']);
    % sgt.FontSize = 35;
end