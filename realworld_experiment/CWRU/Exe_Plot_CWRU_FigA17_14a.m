clc;clear;close all;
loaddatacwru;
if 0
    tau_1 = 1; %the step between start position 1~fix(tau_2/2)
    tau_2 = 3; %the delay lag 2~20???
    d = 11; %embedding dimension: need a range???
    d_p = 4; %length of permutation 3~fix(num_sample/10)
    parametersval = [tau_1,tau_2,d,d_p];
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
    save test4_cwru_100225.mat
else
    load test4_cwru_100225.mat
end
%"Dissipative Standard Map,DSM","Shaw Van Derpol Oscillator, SVDO", "Rossler Attractor,RA",
% "Henon Area Preserving Quadratic Map, HAPQM","Simplest Driven Chaotic Flow, SDCF",
% "Random Noise, RN","White Gaussian Noise, WGN","Random Walk, RW"

%% plot graph
nrows = 2;ncol = 2;
posimat = figposi(nrows,ncol);linewidth_t = 1;
countf = 13;countdrawout = 14;
fonttxt_size = 10;fig_width = 17.6;fig_height = fig_width*0.5;
drawout = 1;
clc;
num_groups = length(unique(cwrulabel));
titlestr_feature = {'d_{cos}U','d_{cos)V','SV_1','SV_2'};
legendstr1={"Bear Faults","Inner Faults","Normal",'Outer Faults'};
legendstr = {"BF_D","BF_F","IF_D","IF_F","N_D","N_F","OF_D","OF_F"};
legendord = {"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"};

markerset = {"d","s","*","o"};
colorset = hsv(num_groups);
if 1
    figure
     plotboard = tiledlayout(8,1);
     plotboard.TileSpacing = 'compact';
    for i_signal = 1:num_groups

        temp_X_mts = X_mts_cell{i_signal*num_samples};
        for i_ord = 1:2
            nexttile;


            pic_signal = plot(temp_X_mts(i_ord,:));
            % tbl_temp_X_mts = array2table(temp_X_mts',"VariableNames",["acc_d","acc_f"]);

            % pic_signal = stackedplot(tbl_temp_X_mts);
            pic_signal.Marker = '.';
            pic_signal.MarkerSize = 2;
            pic_signal.MarkerEdgeColor = colorset(i_signal,:);
            pic_signal.LineStyle = "-";
            pic_signal.LineWidth = linewidth_t/2;
            pic_signal.Color = colorset(i_signal,:);
            tickflag = i_signal==num_groups && i_ord==2;
            if ~tickflag
                set(gca,'XTick',[],'YTick',[]);
            else
                set(gca,'YTick',[]);
            end
            if 1
                ylabel(legendstr{(i_signal-1)*2+i_ord},'FontName','Times New Roman','FontSize',fonttxt_size,'FontWeight','bold')
            end
        end
    end
end

fig_description = 'CWRU_signals';
lowerfont=0;
plotandprint;
num_features = length(titlestr_feature);

if 0
    figure
    tiledlayout(4,1)
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
            wfp(i_cate).LineWidth = 1;
            wfp(i_cate).Color = colorset(i_cate,:);
            % set(gca,"YLim",[4.5,7])
            if i_feature <4
                set(gca,'XTick',[]);
            end

        end

        if i_feature == 1
            lgd = legend(legendstr);
            lgd.NumColumns = 5;
            lgd.Location = 'best';
            lgd.Box = 'off';
            lgd.Color = 'none';
            lgd.FontName = 'Times New Roman';

        end
    end
end
emb_pt_mts = tsne(pe_mts');
figure
gs_chaosnoise = gscatter(emb_pt_mts(:,1),emb_pt_mts(:,2),cwrulabel);
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
lgd = legend(legendstr1);
lgd.NumColumns = 1;
lgd.Location = 'best';
lgd.Box = 'off';
lgd.Color = 'none';
lgd.FontName = 'Times New Roman';
fig_description = 'CWRU_signals';
drawout = 1;
lowerfont=0;
plotandprint;

Traindata = [pe_mts',cwrulabel];