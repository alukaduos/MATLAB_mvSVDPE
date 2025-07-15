clc;clear;close all;%RN 形态+SVDPE
if 1
    n = 1e4;%length
    m = 4;%dimension
    % X_mts = normrnd(0,1,[m,n]);
    % X_mts = repmat(1:n,[m,1]);

    % plot fig parameters
    countdrawout = 1;
    nrows = 2;ncol = 2;
    posimat = figposi(nrows,ncol);linewidth_t = 1;
    countf = 0;countdrawout = 1;
    

    %SVD Entropy parameters
    tau_1 = 1; %the step between start position 1~fix(tau_2/2)
    tau_2 = 3; %the delay lag 2~20???
    d = 7; %embedding dimension: need a range???
    d_p = 5; %length of permutation 3~fix(num_sample/10)
    parametersval = [tau_1,tau_2,d,d_p];

    %Simulation
    simutimes = 100;
    pe_mts = zeros(m+2,simutimes);

    for i = 1:simutimes
        % X_mts = normrnd(0,1/3,[m,n])*rand(1,1)*5;
        X_mts = 1*rand([m,n])*rand(1,1)*5;
        % X_mts(end,:) = 0*X_mts(1,:)+0.02*rand([1,n]);.*repmat(randi(5,[1,n]),[m,1]);
        pe_mts(:,i) = calsvdpe_mts(X_mts,parametersval);
    end
    save hesrn.mat
else
    load hesrn.mat
end
%% plot graph

titlestr_feature = {'PE_{disU}','PE_{disV}','PE_{SV}^1','PE_{SV}^2','PE_{SV}^3',...
    'PE_{SV}^4','PE_{SV}^5','PE_{SV}^6','PE_{SV}^7'};
drawout = 1;

viewindx = 1:15;
figure
tiledlayout(3,3,"TileSpacing","compact")
for i = 1:9
    nexttile;
    if i ==1
        xview = plot(X_mts(1,:),X_mts(2,:));
        xview.Marker = "o";
        xview.MarkerSize = 1;
        xview.LineStyle = "none";
        xview.LineWidth = 0.1;
        xview.Color = 'r';
        % set(gca,'XLim',[0,5],'YLim',[0,5]);
    else
        xview = plot(X_mts(1,viewindx+randi(9000)),X_mts(2,viewindx+randi(9000)));
        xview.Marker = "o";
        xview.MarkerSize = 4;
        xview.LineStyle = ":";
        xview.LineWidth = 1;
        xview.Color = 'r';
        %set(gca,'XLim',[0,5],'YLim',[0,5]);
    end
    set(gca,'XTick',[],'YTick',[])
end
fig_width = 8;fig_height = fig_width;
fig_description = '0_Random_Noise';
fonttxt_size = 8;lowerfont=0;
plotandprint;

countf = countf-1;
countdrawout = countdrawout -1;
% countf = countf+1;
drawout = 1;
singleext = 3;
    bothext = singleext*2;
    viocolorset = hot(6+bothext);
figure
%xgroupdata = {'PE_{disU}^{cos}','PE_{disV}^{cos}','PE_{SV}^1','PE_{SV}^2','PE_{SV}^3', 'PE_{SV}^4'};
bc1 = boxchart(pe_mts','BoxEdgeColor','r','BoxFaceColor','r');
bc1.MarkerSize = 2;
bc1.MarkerColor = 'r';%viocolorset(1,:);
bc1.WhiskerLineColor = bc1.MarkerColor;
ylim([5.5,7])

% xlabel(titlestr_feature);
set(gca,"XTickLabel",titlestr_feature,'Box','on')
ax = gca; 
ax.YTick = 5.5:0.1:7;
grid on
fig_description = '1_SVD_Entropy_for_Random_Noise';
fig_width = 4.1;fig_height = fig_width/0.618;fonttxt_size = 8;lowerfont=1;

if 0
    xgroupdata = {'PE_{disU}^{cos}','PE_{disV}^{cos}','PE_{SV}^1','PE_{SV}^2','PE_{SV}^3', 'PE_{SV}^4'};
    
    hold on
    vio1 = violinplot(pe_mts');
    for i =1:6
        vio1(i).EdgeColor = viocolorset(i+singleext-2,:);
        vio1(i).FaceColor = viocolorset(i+singleext-2,:);
         
    end
    hold off
end


plotandprint;
%%
% countf = countf+1;

% waterfall(pe_mts)
%
% function pe_mts = calsvdpe_mts(X_mts,parametersval)
% tau_1 = parametersval(1);tau_2 = parametersval(2);d = parametersval(3); d_p = parametersval(4);
% [m,n] = size(X_mts);
% k_max = fix((n - (d-1)*tau_2)/tau_1);
% Y = zeros(m,k_max,d);%3D-rearrangements of X_mts, dimension*k_max*emd_dim
% for k = 1:k_max
%     for j = 0:d-1
%         Y(:,k,j+1) = X_mts(:,k*tau_1+j*tau_2);
%     end
% end
% Y_p = permute(Y,[1,3,2]); %dimension*emd_dim*k_max
% s_v = squeeze(pagesvd(Y_p)); %singular val of Y_p
% [U_v,~,V_v] = pagesvd(Y_p); %orthogonal U and V
% U_dim = size(U_v,1); V_dim = size(V_v,1);dis_cos_u = zeros(1,k_max); dis_cos_v = dis_cos_u;
% for k = 1:k_max
%     temp_u = U_v(:,:,k);
%     mean_u = mean(temp_u,2);
%     mean_u_eye = sum(eye(U_dim),2)/U_dim;
%     temp_v = V_v(:,:,k);
%     mean_v = mean(temp_v,2);
%     mean_v_eye = sum(eye(V_dim),2)/V_dim;
%     dis_cos_u(k) = dot(mean_u, mean_u_eye)/norm(mean_u)/norm(mean_u_eye);
%     dis_cos_v(k) = dot(mean_v, mean_v_eye)/norm(mean_v)/norm(mean_v_eye);
%
% end
% S_v = [dis_cos_u;dis_cos_v;s_v];
% i_max = size(S_v,2)-d_p+1;%
% m_p = min(m,d);
% perm = zeros(m_p+2,d_p,i_max);rank_perm = perm;
% for i = 1:i_max
%     perm(:,:,i) = S_v(:,i:i+d_p-1);
%     rank_perm(:,:,i) = tiedrank(perm(:,:,i)')';
% end
% pe_mts = calpe_permat(rank_perm);
% end
% function pe_val = calpe_permat(x)%each row a permutation, rows=(dis_U;dis_V;S_1;S_2;......); each layer an embedded mat
% num_feature = size(x,1);
% pe_val = zeros(num_feature,1);
% for i = 1:num_feature
%     temp_x = squeeze(x(i,:,:));
%     pe_val(i) = calpe_permarray(temp_x);
% end
% end
% function pe_val = calpe_permarray(x)
% [len_perm,num_sample] = size(x);
% [~,~,idx_perm] = unique(x','rows');
% idx_counts = accumarray(idx_perm,1);
% perm_ratio = idx_counts/size(x,2);
% pe_val = -sum(log2(perm_ratio).*perm_ratio);
% end
