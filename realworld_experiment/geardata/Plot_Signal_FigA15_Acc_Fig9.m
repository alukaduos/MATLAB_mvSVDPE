clc;clear; close all;
load geardata_5.mat;
legendstr={"Chipped","Health","Miss",'Root','Surface'};
legendord = {"(a)","(b)","(c)","(d)","(e)"};
% markerset = {"^","+","s","d","o"};
colorset = hsv(num_groups);

nrows = 1;ncol = 2;
posimat = figposi(nrows,ncol);linewidth_t = 1;
countf = 8;countdrawout = countf+1;
fonttxt_size = 8;fig_width = 14;fig_height = fig_width*0.618*0.5;
drawout = 1;
colname = {'MV','PGV_x','PGV_y','PGV_z'};
lowerfont = 0;
% show the signals
if 1
    figure
    dimsignal = length(colname);
    plotboard = tiledlayout(dimsignal,num_groups);
    plotboard.TileSpacing = 'compact';
    for i_ord = 1:dimsignal


        for i_signal = 1:num_groups
            temp_X_mts = X_mts_cell{i_signal*num_samples};
            temp_X_mts = temp_X_mts(:,1:10:end);
            nexttile;
            pic_signal = plot(1:10:8000,temp_X_mts(i_ord,:));
            pic_signal.Marker = 'o';
            pic_signal.MarkerSize = 1;
            pic_signal.MarkerEdgeColor = colorset(i_signal,:);
            pic_signal.LineStyle = "-";
            pic_signal.LineWidth = linewidth_t/2;
            pic_signal.Color = colorset(i_signal,:);
            xtickflag = i_ord==dimsignal;
            % ytickflag = i_signal==1;
            ylabelflag = i_signal==1;
            if ~xtickflag
                set(gca,'XTick',[],'YTick',[]);
            else
                set(gca,'YTick',[]);
            end
            if ylabelflag
                ylabel(colname{i_ord})
            end
            set(gca,"FontSize",fonttxt_size-3,"FontName",'Times New Roman','FontWeight','bold')
            if i_ord == 1
                title(legendord{i_signal},"FontSize",fonttxt_size-3,"FontName",'Times New Roman','FontWeight','bold')
            end

        end
    end

end

fig_description = 'Gear_dataset';
fonttxt_size = fonttxt_size-3;
plotandprint;





load traindata_gear_5.mat
drawout = 1;fig_width = 14; fig_height = fig_width*0.618*0.8;
fonttxt_size = 8;
if 1
    figure
    tiledlayout(nn1,nn2,"TileSpacing","compact")
    colorbarflag = 0;
    for i_1d = 1:nn1
        for i_2d = 1:nn2
            nexttile
            colorbarflag = i_2d==nn2;
            temphp = heatmap(squeeze(trainaccuracy(i_1d,i_2d,:,:)),'Colormap',jet(15));
            % temphp.Colormap = Turbo;
            temphp.ColorLimits = [0.71, 0.91];
            if colorbarflag
                temphp.ColorbarVisible = 'on';
            else
                temphp.ColorbarVisible = 'off';
            end
            temphp.GridVisible = 'off';
            ax = gca;
            ax.XDisplayLabels = nan(size(ax.XDisplayData));
            ax.YDisplayLabels = nan(size(ax.YDisplayData));

            if i_2d == 1
                ax = gca;
                ax.YDisplayLabels = num2cell(5:20);
                ax.YLabel = 'd';
            end
            if i_1d == nn1
                ax = gca;
                ax.XDisplayLabels = num2cell(3:10);
                ax.XLabel = 'd_p';
            end
            temptitle = strcat('\tau_1=',num2str(tau_1_set(i_1d)),', \tau_2=',num2str(tau_2_set(i_2d)));
            set(gca,"FontSize",fonttxt_size-4,"FontName",'Times New Roman')
            title(temptitle);
            % titles.FontSize=fonttxt_size-3;
            % titles.FontName ='Times New Roman';

            % temphp.XDisplayData = [];
            % temphp.YDisplayData = [];
        end
    end
    fig_description = 'TrainAccuracy';
fonttxt_size = fonttxt_size-4;
    
    plotandprint;
end