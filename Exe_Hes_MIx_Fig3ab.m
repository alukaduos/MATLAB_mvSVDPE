clc;clear;close all; % Mixed random signals, coorelation level
if 0
    n = 10e3;%length
    m = 2;%dimension
    


    tau_1 = 1; %the step between start position 1~fix(tau_2/2)
    tau_2 = 3; %the delay lag 2~20???
    d = 7; %embedding dimension: need a range???
    d_p = 5; %length of permutation 3~fix(num_sample/10)
    parametersval = [tau_1,tau_2,d,d_p];

    % chaosdataprocess;
    stochosdataprocess;




   
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
                0.5+((i_stochastic-1)/3+j_stochastic/simutimes/3)/2
            end
        end
        pe_mts = [pe_mts,pe_mts_stochastic];
    end
    save hesmix.mat
else
    load hesmix.mat
    
end
pe_mts_folded = reshape(pe_mts,4,100,[]);
    pe_mts_f_mean = squeeze(mean(pe_mts_folded,2));
    pe_mts_f_std  = squeeze(std(pe_mts_folded,0,2));
%"Dissipative Standard Map,DSM","Shaw Van Derpol Oscillator, SVDO", "Rossler Attractor,RA",
% "Henon Area Preserving Quadratic Map, HAPQM","Simplest Driven Chaotic Flow, SDCF",
% "Random Noise, RN","White Gaussian Noise, WGN","Random Walk, RW"

%% plot graph
lowerfont=0;
    legendstr={"RN","WGN","RmW","Rf","Wf","-0.8","-0.3","0.3","0.8"};
    legendstra=["RN","WGN","RmW","Rf","Wf","-0.8","-0.3","0.3","0.8"];

    ylabelstr = {'PE_{disU}','PE_{disV}','PE_{SV}^1','PE_{SV}^2'};
    markerset = {"+","*","x",".","|",'o','s','v','h'};

nrows = 2;ncol = 2;
posimat = figposi(nrows,ncol);linewidth_t = 1;
%已绘制图片数 countf，下一张保存图片编号 countdrawout
countf = 3;countdrawout = 4;
fonttxt_size = 12;fig_width = 14;fig_height = fig_width*0.618;
% show the data
legendord={"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"};
figure
drawout = 1;
tiledlayout(2,2,'TileSpacing','tight')
colorset = distinguishable_colors(9);
for i = 1:4
    nexttile;
    tempsvdpe = squeeze(pe_mts_folded(i,:,:));
    xgroupdata = categorical(legendstra);
    temppic = boxchart(tempsvdpe);
        temppic.LineWidth = 0.2;
        temppic.MarkerSize = 1;
    box on
    grid on
    ylabel(ylabelstr{i});
    set(gca,"FontSize",fonttxt_size-3,"FontName",'Times New Roman','FontWeight','bold')
    if i>=3
    set(gca,'XTickLabel',legendstr)
    else
        set(gca,'XTickLabel',[])
    end

    if 0
    temppic = errorbar(pe_mts_f_mean(i,1:9),pe_mts_f_std(i,1:9));
    temppic.Color = 0*colorset(i,:);
    temppic.LineWidth = 1;
    ylabel(ylabelstr{i});
    set(gca,"FontSize",fonttxt_size-3,"FontName",'Times New Roman','FontWeight','bold')
    if i>=4
    set(gca,'XTickLabel',legendstr)
    else
        set(gca,'XTickLabel',[])
    end
    end


end
if 0
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
    pic_signal.LineStyle = "-";
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
end
fig_description = '1_Dependent_Signals';
lowerfont=1;
plotandprint;

lowerfont = 0;
emb_pt_mts = tsne(pe_mts(:,1:900)');
mtslable = mtslable(1:900);
colorset = distinguishable_colors(9);

%%
countf = 3;countdrawout = 4;
drawout = 1;
figure
gs_chaosnoise = gscatter(emb_pt_mts(:,1),emb_pt_mts(:,2),mtslable,colorcube(length(unique(mtslable))));
num_catelogs = length(gs_chaosnoise); colorset = distinguishable_colors(num_catelogs); 
%colorset(7,:) = colorset(8,:)-0.18;
% colorset(8,:) = colorset(6,:);
%colorset(6,:) = colorset(9,:)-0.18;
for i_cate = 1:num_catelogs
    gs_chaosnoise(i_cate).Color = colorset(i_cate,:);
    gs_chaosnoise(i_cate).LineWidth = 0.5;
    gs_chaosnoise(i_cate).Marker = markerset{i_cate};
    gs_chaosnoise(i_cate).MarkerSize = 3;
    if i_cate >=8
        gs_chaosnoise(i_cate).MarkerFaceColor= colorset(i_cate,:);
    end
end
set(gca,'XTickLabel',[],'YTickLabel',[])
lgd = legend(legendstr);
lgd.NumColumns = 1;
lgd.Location = 'best';
lgd.Box = 'off';
lgd.Color = 'none';
lgd.FontName = 'Times New Roman';
%lgd.FontSize = 8;
lgd.Location = 'eastoutside';
fig_description = '2_Discremination_Dependent_Signals';
fig_width = 14;fig_height = fig_width*0.618; 
plotandprint;

%Traindata = [pe_mts',labels];



