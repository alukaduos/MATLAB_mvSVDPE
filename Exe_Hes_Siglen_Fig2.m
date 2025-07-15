clc;clear;close all;% WGN different length
if 1
    n_set = [1e3,2e3,5e3,1e4];%length

    m = 4;%dimension
  

    % plot fig parameters
    nrows = 4;ncol = 4;
    posimat = figposi(nrows,ncol);linewidth_t = 1;



    %SVD Entropy parameters
    tau_1 = 1; %the step between start position 1~fix(tau_2/2)
    tau_2 = 3; %the delay lag 2~20???
    d = 7; %embedding dimension: need a range???
    d_p = 5; %length of permutation 3~fix(num_sample/10)
    parametersval = [tau_1,tau_2,d,d_p];

    %Simulation

    simutimes = 100;
    tbltest = table('Size',[length(n_set)*(m+2)*simutimes,3],'VariableTypes',["categorical","double","double"],'VariableNames',{'fea','mts_length','svden'});
    pe_mts = zeros(m+2,simutimes);
    pe_mts_cell = cell(length(n_set),1);
    countsamples = 1;
    for idx_n = 1:length(n_set)
        n = n_set(idx_n);
        for i = 1:simutimes
            X_mts = normrnd(0,1,[m,n]);
            % X_mts = rand([m,n]);
            % X_mts(end,:) = 0*X_mts(1,:)+0.02*rand([1,n]);
            pe_mts(:,i) = calsvdpe_mts(X_mts,parametersval);
            tbltest.mts_length(countsamples:countsamples+m+1) = n_set(idx_n);
            tbltest.fea(countsamples:countsamples+m+1) = categorical({'PE_{disU}','PE_{disV}','PE_{SV}^1',...
                'PE_{SV}^2','PE_{SV}^3', 'PE_{SV}^4'}');
            tbltest.svden(countsamples:countsamples+m+1) = pe_mts(:,i);
            countsamples = countsamples+m+2;
        end
    end
    save hessiglen.mat
else
    load hessiglen.mat
end

%% plot graph
countf = 2;countdrawout = 3;
fonttxt_size = 12;fig_width = 17.6;fig_height = fig_width*0.382;
titlestr_feature = {'PE_{disU}','PE_{disV}','PE_{SV}^1','PE_{SV}^2','PE_{SV}^3',...
    'PE_{SV}^4','PE_{SV}^5','PE_{SV}^6'};
feaorder = {'PE_{disU}','PE_{disV}','PE_{SV}^1','PE_{SV}^2','PE_{SV}^3','PE_{SV}^4'};
tbltest.fea = categorical(tbltest.fea,feaorder);
legendord = {"1000","2000","5000","10000"};
drawout = 1;
figure
boxc1 = boxchart(tbltest.fea,tbltest.svden,'GroupByColor',tbltest.mts_length);
for i_boxc = 1:4
    boxc1(i_boxc).MarkerSize = 1;
    boxc1(i_boxc).LineWidth = 0.5;
end

lgd = legend(legendord);
% lgd = legend(legendstr(1:5));
lgd.NumColumns = 2;
lgd.Location = 'best';
lgd.Box = 'off';
lgd.Color = 'none';
lgd.FontName = 'Times New Roman';
lgd.FontSize = 8;
lgd.FontWeight = "normal";
ylim([5.8,7])
grid on
set(gca,'Box','on')

fig_description = 'SVD_Entropy_for_Different_MTS_Length';
lowerfont = 0;
plotandprint;
%%
