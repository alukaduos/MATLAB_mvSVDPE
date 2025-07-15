clc;clear;close all;%rossler chaos vs limit orbit
recal = 0;% first time exe let recal = 1: 1 = to compute chaos data, 0 = no compute
if 0 % first time exe let this be 1: take long time to compute
    %Time Series parameters
    series_length = 10000; time_step = 1e-2; bifurcation = 0.37:0.0001:0.43;
    %SVD Entropy parameters
    tau_1 = 1; %the step between start position 1~fix(tau_2/2)
    tau_2 = 2; %the delay lag 2~20???
    d = 7; %embedding dimension: need a range???
    d_p = 5; %length of permutation 3~fix(num_sample/10)
    parametersval = [tau_1,tau_2,d,d_p];

    %Initial data
    if recal
        num_signals = length(bifurcation);
        X_mts_cell = cell(num_signals,1);
        xweizhi = X_mts_cell;
        pe_mts = zeros(5,num_signals);
        ly_mts = zeros(1,num_signals);
        for i_signal = 1:num_signals
            X_mts_cell{i_signal} = sim_rosslerattractor_mts_b(series_length,time_step,bifurcation(i_signal))';
            i_signal/num_signals
        end
        save chaoslimitsimdata.mat
    else
        load chaoslimitsimdatasimdata.mat
    end
    %ly_mts = zeros(1,num_signals);
    extractgap = 3;
    extractindx = 1:extractgap:30000*extractgap+1;
    for i_signal = 1:num_signals
        temp_lya = zeros(3,1);
        temp_Xmts = X_mts_cell{i_signal}(1,:);
        temp_Xmts = temp_Xmts(:,extractindx);
        index_x = abs(diff(temp_Xmts))<=0.003;
        xweizhi{i_signal} = unique(temp_Xmts(index_x));
        pe_mts(:,i_signal) = calsvdpe_mts(X_mts_cell{i_signal}(:,extractindx),parametersval);
        % for j = 1:3
        % [~,~,temp_lya(j)] = chaostest(X_mts_cell{i_signal}(j,extractindx));
        % end
        % ly_mts = max(temp_lya);

        i_signal/num_signals
    end
    spoints = [0.3906, 0.3908, 0.3934, 0.3964, 0.3997, 0.3999, 0.4032, 0.4088,0.4183,0.4245] ;
    save chaoslimitsimdata.mat
else
    load chaoslimitsimdata.mat
end
chaoticindx = find(ly_mts<=0);
chaotica = bifurcation(chaoticindx);
showca = chaotica(chaotica>=0.3906);

%% plot
nrows = 2;ncol = 2;
posimat = figposi(nrows,ncol);linewidth_t = 1;
countf = 7;countdrawout = 8;
titlestr_feature = {'PE_{disU}','PE_{disV}','PE_{SV}^1','PE_{SV}^2','PE_{SV}^3'};
legendord={"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"};

fonttxt_size = 8;fig_width = 17.6;fig_height = fig_width*0.5;
drawout = 1;
figure
tiledlayout(3,2,'TileSpacing','compact')
for i_pic = 1:6
    nexttile;
    hold on
    if i_pic ==1
        for i_signal = 1:num_signals
            plot(bifurcation(i_signal)*ones(length((xweizhi{i_signal})),1), (xweizhi{i_signal}), 'b.', 'markersize',1)
        end
        plotverticallines(spoints,[-4,6]);
        hold off
        axis([0.37, 0.43,1,5.5])
    else
        plot(bifurcation,pe_mts(i_pic-1,:),'bo--','markersize',2,'LineWidth',0.2)
        hold on
        plotverticallines(spoints,[min(pe_mts(i_pic-1,:)),max(pe_mts(i_pic-1,:))]);
        hold off
        ylabel(titlestr_feature{i_pic-1},'FontName','Times New Roman','FontSize',fonttxt_size,'FontWeight','bold');
        axis([0.37, 0.43,min(pe_mts(i_pic-1,:))*0.999,max(pe_mts(i_pic-1,:))*1.001])
    end
    set(gca,"FontSize",fonttxt_size,"FontName",'Times New Roman','FontWeight','bold')
    title(legendord{i_pic},'FontName','Times New Roman','FontSize',fonttxt_size,'FontWeight','bold');
    if i_pic<5
        xticklabels('')
    else
        xticks(0.37:0.01:0.43)
    end
    %yticklabels('')
    box on
    grid on
end
fig_description = 's_Rossler_chaosorbit';
lowerfont = 0;
plotandprint;



function plotverticallines(spoints, data)
n = length(spoints);
maxp = max(data)*1.001;minp = min(data)*0.997;
for i=1:n
    if ismember(i,[1,9,10])
        plot([spoints(i), spoints(i)], [maxp,minp],'r:','linewidth',0.3)
    else
        plot([spoints(i), spoints(i)], [maxp,minp],'r:','linewidth',0.3)
    end
end
end




function y = sim_rosslerattractor_mts_b(series_length,time_step,birf)
global a
a = birf;
n = fix(1.5*series_length);
gapscale = 10;
gap = 1/time_step/gapscale;
tspan = 0:time_step:n;
initial_sim_y = [1 0 1]';
[~,sim_y] = ode45(@rossler_z,tspan,initial_sim_y);
clear a;
% plot(sim_y(:,1),sim_y(:,2),'o','MarkerSize',1)

sampleindx = 1:gap:length(sim_y(:,1));
sim_y = sim_y(sampleindx,:);
sim_y = zscore(sim_y);

y = sim_y(end-series_length*gapscale+1:end,:);

    function dy = rossler_z(t,y)
        % global a
        b = 2; c = 4;
        dy = zeros(3,1);
        dy(1) = -y(2) - y(3);
        dy(2) = y(1) + a*y(2);
        dy(3) = b + y(3)*(y(1) - c);
    end
end


