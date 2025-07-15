clc;clear;close all;
if 1
    loadgeardata;
    % load geardata_5.mat
    tau_1 = 1; %the step between start position 1~fix(tau_2/2)
    tau_2 = 11; %the delay lag 2~20???
    d = 5; %embedding dimension: need a range???
    d_p = 4; %length of permutation 3~fix(num_sample/10)
    parametersval = [tau_1,tau_2,d,d_p];
    pe_mts_cell = cell(num_totalsamples,1);
    oripe_mts_cell  = pe_mts_cell;
    dim_pe_mts = zeros(num_totalsamples,1);



    for i = 1:num_totalsamples
        % X_mts = normrnd(0,2,[m,n]);
        X_mts = X_mts_cell{i};
        % X_mts(end,:) = 0*X_mts(1,:)+0.02*rand([1,n]);
        pe_mts_cell{i} = calsvdpe_mts(X_mts,parametersval);
        oripe_mts_cell{i} = calpe_mts(X_mts,parametersval);
        dim_pe_mts(i) = length(pe_mts_cell{i});
        i/num_totalsamples/2
    end
    dim_pe_vec = max(dim_pe_mts);
    pe_mts = zeros(dim_pe_vec,num_totalsamples);
    oripe_mts = zeros(dim_pe_vec-2,num_totalsamples);
    for i = 1:num_totalsamples
        if dim_pe_mts(i)<dim_pe_vec
            pe_mts(:,i) = padarray(pe_mts_cell{i},dim_pe_vec-dim_pe_mts(i),0,'post');
            oripe_mts(:,i) = padarray(oripe_mts_cell{i},dim_pe_vec-dim_pe_mts(i),0,'post');
        else
            pe_mts(:,i) = pe_mts_cell{i};
            oripe_mts(:,i) = oripe_mts_cell{i};
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
    save test4_gear_110225.mat
else
    load test4_gear_110225.mat
end

%% plot graph
clc;close all;
num_groups = length(unique(gearlabel));
titlestr_feature = {'d_{cos}U','d_{cos)V','SingularValue_1','SingularValue_2','SingularValue_3','SingularValue_4',...
    'SingularValue_5','SingularValue_6','SingularValue_7'};
legendstr={"Chipped","Health","Miss",'Root','Surface'};
markerset = {"^","+","s","d","o"};
colorset = jet(num_groups);
if 0
    figure
    tiledlayout(2,5)
    for i_signal = 1:num_groups

        nexttile;
        temp_X_mts = X_mts_cell{i_signal*num_samples};
        pic_signal = plot(temp_X_mts(1,:),temp_X_mts(2,:));
        pic_signal.Marker = '.';
        pic_signal.MarkerEdgeColor = colorset(i_signal,:);
        pic_signal.LineStyle = "none";
        pic_signal.Color = colorset(i_signal,:);
        set(gca,'XTick',[],'YTick',[])
        title(legendstr{i_signal})
    end
end
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
if 0
figure
gs_chaosnoise = gscatter(emb_pt_mts(:,1),emb_pt_mts(:,2),gearlabel);
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
lgd.NumColumns = 1;
lgd.Location = 'best';
lgd.Box = 'off';
lgd.Color = 'none';
lgd.FontName = 'Times New Roman';
title(num2str(selectcol))
end
Traindata = [pe_mts',gearlabel];
Traindata_ori = [oripe_mts',gearlabel];

% Then use Classifier APP to train and test