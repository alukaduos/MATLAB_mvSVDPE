clc;clear;close all; % rossler+noise+enlarge
if 0
    n = 1e3;%length
    m = 6;%dimension
   


    tau_1 = 1; %the step between start position 1~fix(tau_2/2)
    tau_2 = 3; %the delay lag 2~20???
    d = 7; %embedding dimension: need a range???
    d_p = 5; %length of permutation 3~fix(num_sample/10)
    parametersval = [tau_1,tau_2,d,d_p];

    rossleraddnoise;

   
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
    save chaos_noise.mat
else
    load chaos_noise.mat
end
%"Dissipative Standard Map,DSM","Shaw Van Derpol Oscillator, SVDO", "Rossler Attractor,RA",
% "Henon Area Preserving Quadratic Map, HAPQM","Simplest Driven Chaotic Flow, SDCF",
% "Random Noise, RN","White Gaussian Noise, WGN","Random Walk, RW"

%% plot graph
nrows = 2;ncol = 2;lowerfont=0;
posimat = figposi(nrows,ncol);linewidth_t = 1;
%已绘制图片数 countf，下一张保存图片编号 countdrawout
countf = 5;countdrawout = 6;
fonttxt_size = 8;fig_width = 8;fig_height = fig_width;
titlestr_feature = {'d_{cos}U','d_{cos)V','SingularValue_1','SingularValue_2','SingularValue_3'};
legendstr={"original","0.03","0.05","0.10","0.15","0.20","0.30","0.50",'wgn'};
legendord={"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"};

markerset = {".","o","s","d","^",'v','<','>',"p"};
colorset = jet(num_groups);

drawout = 1;
figure
tiledlayout(3,3)
for i_signal = 1:num_groups

    nexttile;
    temp_X_mts = X_mts_cell{i_signal*100};
    pic_signal = plot3(temp_X_mts(1,:),temp_X_mts(2,:),temp_X_mts(3,:));
    pic_signal.Marker = 'd';
    pic_signal.MarkerSize = 3;
    pic_signal.MarkerEdgeColor = colorset(i_signal,:);
    pic_signal.LineStyle = ":";
    pic_signal.Color = colorset(i_signal,:);
    pic_signal.LineWidth = 0.02;
    set(gca,'XTick',[],'YTick',[],'ZTick',[])
    set(gca,"FontSize",fonttxt_size,"FontName",'Times New Roman','FontWeight','bold')
    title(legendord{i_signal},'FontName','Times New Roman','FontSize',fonttxt_size,'FontWeight','bold')
    set(gca,"Box","on")
    grid on
end
fig_description = 'Rossler_with_Noise';
plotandprint;

num_features = length(titlestr_feature);

% fig_width = 14; fig_height = fig_width*0.618;



figure
tiledlayout(3,2,"TileSpacing","tight")
for i_feature = 1:num_features+1
    if i_feature == 6
        lgd = legend(legendstr(1:9));
        lgd.NumColumns = 3;
        lgd.Location = 'best';
        lgd.Box = 'off';
        lgd.Color = 'white';
        lgd.FontName = 'Times New Roman';
        lgd.FontSize = 12;
        lgd.FontWeight = "normal";
        lgd.Layout.Tile = 6;  
        lgd.IconColumnWidth = 18;

        break

    else
        nexttile;
    end

    temp_featureval = pe_mts(i_feature,:)';
    temp_featureval = reshape(temp_featureval,num_samples,[]);
    wfp = plot(temp_featureval);
    num_catelogs = length(wfp);

    for i_cate = 1:num_catelogs
        wfp(i_cate).Marker = markerset{i_cate};
        wfp(i_cate).MarkerEdgeColor = colorset(i_cate,:);
        wfp(i_cate).MarkerSize = 2;
        wfp(i_cate).LineWidth = 0.5;
        wfp(i_cate).LineStyle = ":";
        wfp(i_cate).Color = colorset(i_cate,:);
        if i_feature<=2
            set(gca,"YLim",[6.2,6.9])
        else
            set(gca,"YLim",[4.5,6.6])
            set(gca,"YTick",4.5:0.3:6.6)
        end

        set(gca,"FontSize",fonttxt_size,"FontName",'Times New Roman','FontWeight','bold')

        if 0&&(i_feature ==2 || i_feature > 3)
            set(gca,'YTick',[]);
            set(gca,'XTick',[]);
        end

    end

    if i_feature == 15
        lgd = legend(legendstr(1:9));
        lgd.NumColumns = 2;
        lgd.Location = 'southoutside';
        lgd.Box = 'off';
        lgd.Color = 'none';
        lgd.FontName = 'Times New Roman';
        lgd.FontSize = 10;
        lgd.FontWeight = "normal";
    elseif i_feature == 10

        lgd = legend(wfp([6:9]),legendstr(6:9));
        lgd.NumColumns = 5;
        lgd.Location = 'best';
        lgd.Box = 'off';
        lgd.Color = 'none';
        lgd.FontName = 'Times New Roman';
        lgd.FontSize = 8;
        lgd.FontWeight = "normal";
    end
    set(gca,"Box","on")
    grid on
    tempt1 = title(legendord{i_feature},'FontName','Times New Roman','FontSize',fonttxt_size,'FontWeight','bold');
end
fig_description = '1_SVD_Entropy_of_Rossler_with_Noise';
fig_width = 14; fig_height = fig_width*0.618;
lowerfont = 0;
drawout = 1;
plotandprint;


lowerfont = 0;
emb_pt_mts = tsne(pe_mts');

fig_width = 14; fig_height = fig_width*0.618;
drawout = 1;

figure
gs_chaosnoise = gscatter(emb_pt_mts(:,1),emb_pt_mts(:,2),labels);
num_catelogs = length(gs_chaosnoise);
% colorset = jet(num_groups);

% colorset = distinguishable_colors(num_catelogs);
for i_cate = 1:num_catelogs
    gs_chaosnoise(i_cate).Color = colorset(i_cate,:);
    gs_chaosnoise(i_cate).LineWidth = 0.5;
    gs_chaosnoise(i_cate).Marker = markerset{i_cate};
    gs_chaosnoise(i_cate).MarkerSize = 3;
    set(gca,'XTickLabel',[],'YTickLabel',[])
    set(gca,"Box","on")
end

lgd = legend(legendstr);
lgd.NumColumns = 1;
lgd.Location = 'eastoutside';
lgd.Box = 'off';
lgd.Color = 'none';
lgd.FontName = 'Times New Roman';
%lgd.FontSize = 12;
fig_description = 'Discremination_Rossler_Chaos_at_Different_Noise_Levels';
fonttxt_size = 12;fig_width = 14;fig_height = fig_width*0.618;
plotandprint;
Traindata = [pe_mts',labels];


