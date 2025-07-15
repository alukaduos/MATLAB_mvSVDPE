clc;clear;close all;%restrict RW 形态+SVDPE
if 1
    n = 1e4;%length
    m = 4;%dimension
    % X_mts = normrnd(0,1,[m,n]);
    % X_mts = repmat(1:n,[m,1]);

    % plot fig parameters
    nrows = 2;ncol = 2;
    posimat = figposi(nrows,ncol);linewidth_t = 1;
    countf = 1;countdrawout = 2;
    

    %SVD Entropy parameters
    tau_1 = 1; %the step between start position 1~fix(tau_2/2)
    tau_2 = 3; %the delay lag 2~20???
    d = 7; %embedding dimension: need a range???
    d_p = 5; %length of permutation 3~fix(num_sample/10)
    parametersval = [tau_1,tau_2,d,d_p];

    %Simulation
    simutimes = 100;
    pe_mts = zeros(m+2,simutimes);
    sidelength = 20;%side of restricted square
    for i = 1:simutimes
        X_mts = normrnd(0,1,[m,n])/3.*(repmat(rand(1,n)*5,[m,1]));
        X_mts = limitwalk(X_mts,0*ones(m,1),sidelength*ones(m,1));
        % X_mts = cumsum(X_mts,2);
        % X_mts = rand([m,n]);
        % X_mts(end,:) = 0*X_mts(1,:)+0.02*rand([1,n]);
        pe_mts(:,i) = calsvdpe_mts(X_mts,parametersval);
    end
    switch sidelength
        case 20
            save hesrw_20.mat
        case 5
            save hesrw_5.mat
    end
else
    if 1 %1 = side==20, 0 = side==5
        load hesrw_20.mat
    else
        load hesrw_5.mat
    end
end
%% plot graph
titlestr_feature = {'PE_{disU}','PE_{disV}','PE_{SV}^1','PE_{SV}^2','PE_{SV}^3',...
    'PE_{SV}^4','PE_{SV}^5','PE_{SV}^6','PE_{SV}^7'};
drawout =1;
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
        xview.Color = 'b';
        set(gca,'XLim',[0,5],'YLim',[0,5]);
    else
        xview = plot(X_mts(1,viewindx+randi(9000)),X_mts(2,viewindx+randi(9000)));
        xview.Marker = "o";
        xview.MarkerSize = 4;
        xview.LineStyle = ":";
        xview.LineWidth = 1;
        xview.Color = 'b';
        % set(gca,'XLim',[0,5],'YLim',[0,5]);
    end
    set(gca,'XTick',[],'YTick',[])
end
fig_description = '0_Restrict_Random_Walk';
fonttxt_size = 8;fig_width = 8;fig_height = fig_width;lowerfont=0;
plotandprint;
countf = countf-1;
countdrawout = countdrawout -1;

drawout = 1;
figure
bc1 = boxchart(pe_mts','BoxEdgeColor','b','BoxFaceColor','b');
bc1.MarkerColor = [0.1 0.1 1];
bc1.MarkerSize = 2;
bc1.WhiskerLineColor = bc1.MarkerColor;

ylim([5.5,7])

% xlabel(titlestr_feature);
set(gca,"XTickLabel",titlestr_feature,'Box','on')
ax = gca; 
ax.YTick = 5.5:0.1:7;
grid on
fig_description = '1_SVD_Entropy_for_restrict_Random_Walk_20';
fonttxt_size = 8;fig_width = 3;fig_height = fig_width/0.618;lowerfont=1;
plotandprint;
%%


function y = limitwalk(X,lbd,ubd)
[m,n] = size(X);
y = zeros(m,n);
for i = 1:m
    y(i,:) = singlelimitwalk(X(i,:),lbd(i),ubd(i));
end

end
function y = singlelimitwalk(X,lbd,ubd)
n = length(X);
y = zeros(1,n);
for i = 1:n
    if i == 1
        if X(i) > ubd || X(i) < lbd
            X(1) = (lbd+ubd)/2;
            y(1) = X(1);
        else
            y(1) = X(1);
        end
    end
    if i>1
        tempy = y(i-1) + X(i);
        if tempy < lbd
            tempy = lbd + (lbd-tempy);
        elseif tempy > ubd
            tempy = ubd - (tempy-ubd);

        end
        y(i) = tempy;
    end

end
end

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
