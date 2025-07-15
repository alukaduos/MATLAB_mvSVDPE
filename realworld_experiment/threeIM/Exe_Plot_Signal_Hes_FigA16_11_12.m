clc;clear;close all;
%n = 1e3;%length
%m = 6;%dimension
% X_mts = normrnd(0,1,[m,n]);
% X_mts = repmat(1:n,[m,1]);
if 0

    tau_1 = 2; %the step between start position 1~fix(tau_2/2)
    tau_2 = 5; %the delay lag 2~20???
    d = 13; %embedding dimension: need a range???
    d_p = 8; %length of permutation 3~fix(num_sample/10)
    parametersval = [tau_1,tau_2,d,d_p];

    % chaosaddnoiseprocess;
    %rossleraddnoise;
    loadhfdata;
    num_totalsamples = length(X_mts_cell);
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
    pe_mts_ori = pe_mts;
    for i = 1:num_totalsamples
        if dim_pe_mts(i)<dim_pe_vec
            pe_mts(:,i) = padarray(pe_mts_cell{i},dim_pe_vec-dim_pe_mts(i),0,'post');
        else
            pe_mts(:,i) = pe_mts_cell{i};
        end
        i/num_totalsamples/2+.5
    end
    if 0
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
                % 1+(i_stochastic-1)/
            end
        end
        pe_mts = [pe_mts,pe_mts_stochastic];
    end
    save test4_hf_080225.mat
else
    load test4_hf_080225.mat
end
%"Dissipative Standard Map,DSM","Shaw Van Derpol Oscillator, SVDO", "Rossler Attractor,RA",
% "Henon Area Preserving Quadratic Map, HAPQM","Simplest Driven Chaotic Flow, SDCF",
% "Random Noise, RN","White Gaussian Noise, WGN","Random Walk, RW"

%% plot graph
titlestr_feature = {'PE_{disU}','PE_{disV}','PE^1_{SV}','PE^2_{SV}','PE^3_{SV}'};
legendstr={"Health","Faulty"};
ylabelset={"x","y","z"};
markerset = {"o","d"};
num_groups = max(vib_label)+1;
colorset = hsv(num_groups+2);linealpha = 0.2;
colorset_temp = colorset;
colorset(1,:) = colorset_temp(2,:);
colorset(2,:) = colorset_temp(1,:);

nrows = 2;ncol = 2;
posimat = figposi(nrows,ncol);linewidth_t = 1;
%已绘制图片数 countf，下一张保存图片编号 countdrawout
countf = 10;countdrawout = countf+1;
fonttxt_size = 8;fig_width = 17.6;fig_height = fig_width*0.3;
drawout = 1;
legendord = {"(a)","(b)","(c)","(d)","(e)"};
figure
ylimits = [-3 3;-5,5;-4,2]
tdl1 = tiledlayout(3,2,"TileSpacing","compact");
num_ele_ingroup = num_hf/2;
for j = 1:3
for i_signal = 1:num_groups
   
    nexttile;

    temp_X_mts = X_mts_cell{i_signal*num_ele_ingroup};
    pic_signal = plot(temp_X_mts(j,1:1:end));
    pic_signal.Marker = '.';
    pic_signal.MarkerSize = 1;
    pic_signal.MarkerFaceColor = colorset(i_signal,:);
    pic_signal.MarkerEdgeColor = colorset(i_signal,:);
    pic_signal.LineStyle = ":";
    pic_signal.LineWidth = 0.5;
    pic_signal.Color = [colorset(i_signal,:),linealpha];
    xlim([0,5000]);
    ylim(ylimits(j,:))
    if i_signal == 1
        ylabel(ylabelset{j},'FontName','Times New Roman','FontSize',fonttxt_size,'FontWeight','bold')
    %else
        %set(gca,'YTickLabel',[])
    end
    if j ~=3
    set(gca,'XTicklabel',[])
    else
        set(gca,"XTick",0:1000:5000)
    end
    set(gca,"Box","on","FontSize",fonttxt_size,"FontName",'Times New Roman','FontWeight','bold')
 
    if j == 2
        set(gca,'YTick', -4:2:4,"YTickLabel",-4:2:4)
    end
    if i_signal == 2
        set(gca,'YTickLabel',[])
    end
    if j == 1
    temp1 = title(legendord{i_signal},'FontName','Times New Roman','FontSize',fonttxt_size,'FontWeight','bold')
    end
    %set(gca,"Box","on","FontSize",fonttxt_size,"FontName",'Times New Roman','FontWeight','bold')
    grid on

end
end
lowerfont = 0;
fig_description = 'HF_dataset';
fig_width = 17.6;fig_height = fig_width*0.3;
plotandprint;
% another figure
num_features = length(titlestr_feature);
num_samples = size(pe_mts,2)/num_groups;
figure
tiledlayout(5,1,"TileSpacing","compact")
for i_feature = 1:num_features
    nexttile;

    temp_featureval = pe_mts(i_feature,:)';
    temp_featureval = reshape(temp_featureval,num_samples,[]);
    wfp = plot(temp_featureval);
    num_catelogs = length(wfp);
    for i_cate = 1:num_catelogs
        wfp(i_cate).Marker = markerset{i_cate};
        wfp(i_cate).MarkerEdgeColor = colorset(i_cate,:);
        wfp(i_cate).MarkerSize = 2;
        wfp(i_cate).LineStyle = "none";
        wfp(i_cate).LineWidth = 1;
        wfp(i_cate).Color = colorset(i_cate,:);
        
        
        %set(gca,"YLim",[4.5,7])
        if i_feature <5
            set(gca,'XTick',[]);
        end

    end

    set(gca,"TickLabelInterpreter",'latex',"Box","on","FontSize",fonttxt_size,"FontName",'Times New Roman','FontWeight','bold')
    ylabel(titlestr_feature{i_feature},"FontSize",fonttxt_size,"FontName",'Times New Roman','FontWeight','bold')

    if i_feature == 1
        lgd = legend(legendstr);
        lgd.NumColumns = 5;
        lgd.Location = 'best';
        lgd.Box = 'off';
        lgd.Color = 'none';
        lgd.FontName = 'Times New Roman';

    end
        title(legendord{i_feature},'FontName','Times New Roman','FontSize',fonttxt_size,'FontWeight','bold')
        set(gca,"Box","on")
        grid on
end
fig_description = 'svdEntropy_HF_dataset';
drawout = 1;
fig_width = 17.6;fig_height = fig_width*.5;
plotandprint;


emb_pt_mts = tsne(pe_mts');
labels = vib_label;
figure
gs_chaosnoise = gscatter(emb_pt_mts(:,1),emb_pt_mts(:,2),labels);
num_catelogs = length(gs_chaosnoise);
% colorset = jet(num_groups);

% colorset = distinguishable_colors(num_catelogs);
for i_cate = 1:num_catelogs
    gs_chaosnoise(i_cate).Color = colorset(i_cate,:);
    gs_chaosnoise(i_cate).LineWidth = 1;
    gs_chaosnoise(i_cate).Marker = markerset{i_cate};
    gs_chaosnoise(i_cate).MarkerSize = 6;
    set(gca,'XTick',[],'YTick',[])
end
lgd = legend(legendstr);
lgd.NumColumns = 2;
lgd.Location = 'best';
lgd.Box = 'off';
lgd.Color = 'none';
lgd.FontName = 'Times New Roman';
fig_description = 'Discramination_HF_datasets';
drawout = 1;
fig_width = 17.6;fig_height = fig_width*0.3;
plotandprint;

%%
Traindata = [pe_mts',labels];
idx_in_1 = randperm(80,24);
idx_in_2 = randperm(80,24)+80;
Traindata_1 = Traindata([idx_in_1,idx_in_2],:);
Traindata_2 = Traindata;
Traindata_2([idx_in_1,idx_in_2],:)=[];
save test4_hf_3_7.mat

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
