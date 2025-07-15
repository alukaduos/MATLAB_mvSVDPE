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
    loadhfdata;
    num_totalsamples = length(X_mts_cell);
    pe_mts_cell = cell(num_totalsamples,1);

    dim_pe_mts = zeros(num_totalsamples,1);



    for i = 1:num_totalsamples
        X_mts = X_mts_cell{i};
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

