clc;clear;close all;% 8signals, chaos vs stochastic
if 0
n = 3e3;%length
m = 6;%dimension



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


