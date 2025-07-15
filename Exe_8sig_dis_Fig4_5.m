clc;clear;close all;% 8signals, chaos vs stochastic
if 0
n = 3e3;%length
m = 6;%dimension
% X_mts = normrnd(0,1,[m,n]);
% X_mts = repmat(1:n,[m,1]);


tau_1 = 1; %the step between start position 1~fix(tau_2/2)
tau_2 = 3; %the delay lag 2~20???
d = 7; %embedding dimension: need a range???
d_p = 5; %length of permutation 3~fix(num_sample/10)
parametersval = [tau_1,tau_2,d,d_p];

chaosdataprocess;


legendstr={"DSM","SVDO","RA","HAPQM","SDCF","RN","WGN","RWB"};
markerset = {"+","*","x",".","|",'o','s','>'};

% simutimes = 10;
% pe_mts = zeros(m+2,simutimes);
pe_mts_cell = cell(num_totalsamples,1);
dim_pe_mts = zeros(num_totalsamples,1);


for i = 1:num_totalsamples
    % X_mts = normrnd(0,2,[m,n]);
    X_mts = X_mts_cell{i};

    % X_mts(end,:) = 0*X_mts(1,:)+0.02*rand([1,n]);
    pe_mts_cell{i} = calsvdpe_mts(X_mts,parametersval);
    dim_pe_mts(i) = length(pe_mts_cell{i});
    i/num_totalsamples/2
end
dim_pe_vec = max(dim_pe_mts);
pe_mts = zeros(dim_pe_vec,num_totalsamples);
for i = 1:num_totalsamples
    if dim_pe_mts(i)<dim_pe_vec
        pe_mts(:,i) = padarray(pe_mts_cell{i},dim_pe_vec-dim_pe_mts(i),0,'post');
    else
        pe_mts(:,i) = pe_mts_cell{i};
    end
    i/num_totalsamples/2+.5
end
simutimes = 100;
pe_mts_stochastic = zeros(dim_pe_vec,3*simutimes);
labels = [labels;zeros(3*simutimes,1)];
for i_stochastic = 1:3
    for j_stochastic = 1:simutimes
        idx_temp = (i_stochastic-1)*simutimes+j_stochastic;
        labels(num_totalsamples+idx_temp) = 5+i_stochastic;
        switch i_stochastic
            case 1
                X_mts_temp = rand(3,length_sample_each);
            case 2
                X_mts_temp = normrnd(0,1,[3,length_sample_each]);
            case 3
                X_mts_temp = cumsum(normrnd(0,1,[3,length_sample_each]),2);
        end
        pe_mts_stochastic(:,idx_temp) = calsvdpe_mts(X_mts_temp,parametersval);
        0.5+((i_stochastic-1)/3+j_stochastic/simutimes/3)/2
    end
end
pe_mts = [pe_mts,pe_mts_stochastic];
save chaosvsstocha.mat
else
    load chaosvsstocha.mat
end
%"Dissipative Standard Map,DSM","Shaw Van Derpol Oscillator, SVDO", "Rossler Attractor,RA",
% "Henon Area Preserving Quadratic Map, HAPQM","Simplest Driven Chaotic Flow, SDCF",
% "Random Noise, RN","White Gaussian Noise, WGN","Random Walk, RW"

%% plot graph
legendstr={"DSM","SVDO","RA","HM","DCF","RN","WGN","RW"};
nrows = 2;ncol = 2;
posimat = figposi(nrows,ncol);linewidth_t = 1;
%已绘制图片数 countf，下一张保存图片编号 countdrawout
countf = 3;countdrawout = 4;
fonttxt_size = 8;fig_width = 17.6;fig_height = fig_width*0.5;
% show the data
legendord={"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"};
figure
drawout = 1;
tiledlayout(2,4)
colorset = distinguishable_colors(8);
for i = 1:8
    nexttile
    if i<=5
        x_mts_plot = X_mts_cell{i*100};


    elseif i==6
        x_mts_plot = rand(3,length_sample_each);
    elseif i==7
        x_mts_plot = normrnd(0,1,[3,length_sample_each]);
    else
        x_mts_plot = cumsum(normrnd(0,1,[3,length_sample_each]),2);
    end

    if i<=5

        switch dim_mts(i*100)
            case 2
                pic_signal = plot(x_mts_plot(1,:),x_mts_plot(2,:));
            case 3
                pic_signal = plot3(x_mts_plot(1,:),x_mts_plot(2,:),x_mts_plot(3,:));
        end
    else
        pic_signal = plot3(x_mts_plot(1,:),x_mts_plot(2,:),x_mts_plot(3,:));
    end
    pic_signal.Marker = markerset{i};
    pic_signal.MarkerEdgeColor = colorset(i,:);
    pic_signal.MarkerSize = 3;
    pic_signal.LineWidth = 0.02;
    pic_signal.LineStyle = ":";
    
    pic_signal.Color = colorset(i,:);
    if i == 3 || i>=6
        set(gca,'XTick',[],'YTick',[],'ZTick',[])
    else
        set(gca,'XTick',[],'YTick',[])
    end
    set(gca,"Box","on")
    grid on
    title(legendord{i},'FontName','Times New Roman','FontSize',fonttxt_size,'FontWeight','bold');
end
fig_description = '8_Kind_of_Signals';
lowerfont = 0;
plotandprint;


emb_pt_mts = tsne(pe_mts');

drawout = 1;
figure
gs_chaosnoise = gscatter(emb_pt_mts(:,1),emb_pt_mts(:,2),labels,colorcube(length(unique(labels))));
num_catelogs = length(gs_chaosnoise); colorset = distinguishable_colors(num_catelogs);
for i_cate = 1:num_catelogs
    gs_chaosnoise(i_cate).Color = colorset(i_cate,:);
    gs_chaosnoise(i_cate).LineWidth = 0.5;
    gs_chaosnoise(i_cate).Marker = markerset{i_cate};
    gs_chaosnoise(i_cate).MarkerSize = 3;
end
set(gca,'XTickLabel',[],'YTickLabel',[])
lgd = legend(legendstr);
lgd.NumColumns = 1;
lgd.Location = 'best';
lgd.Box = 'off';
lgd.Color = 'none';
lgd.FontName = 'Times New Roman';
lgd.Location = 'eastoutside';
fig_description = 'Discremination_8_Signals';
%fig_width = 14;fig_height = fig_width*0.618;
fonttxt_size = 12;fig_width = 17.6;fig_height = fig_width*0.618;

plotandprint;

%Traindata = [pe_mts',labels];


% figure
% boxchart(pe_mts_cell')
%% functions comments
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
