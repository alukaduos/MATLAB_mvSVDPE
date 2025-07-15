clc;clear;close all;
load result_hf_090425.mat
num_classifer = size(result_hf_090425,1);
x = result_hf_090425.Classifier;
y = [result_hf_090425.ra1,result_hf_090425.rt1,result_hf_090425.ra2, result_hf_090425.rt2];
colorset = colorcube(33); 
markersizes = 4;
%colorset = [0.8500 0.3250 0.0980;1 0 1;0.4660 0.6740 0.1880;0 0 1]
% templine = 50:0.1:100;
countf = 12;countdrawout = 13;
nrows = 2;ncol = 2;
posimat = figposi(nrows,ncol);linewidth_t = 1;
fonttxt_size = 10;fig_width = 17.6;fig_height = fig_width*0.382;
drawout = 1;
figure

barplot = bar(x,y,'FaceColor','flat');
for i = 1:4
        barplot(i).FaceColor = colorset(7*i,:);    
end
ylim([90,100])
lgd = legend("RA for H_{ES}","RT for H_{ES}",...
    'RA for H_{P}',"RT for H_{P}");
lgd.Box = "off";
lgd.Location = "northoutside";
lgd.NumColumns = 4;
lgd.Interpreter = "tex";
% 
% hold on
% svdpeplot = plot(y(:,1),y(:,2));
% %svdpeplot.Color = colorset(1,:);
% svdpeplot.Marker = 'o';
% svdpeplot.LineStyle = 'none';
% svdpeplot.MarkerSize = markersizes;
% svdpeplot.MarkerEdgeColor = colorset(1,:);
% oripeplot = plot(y(:,3),y(:,4));
% oripeplot.Marker = '*';
% oripeplot.LineStyle = 'none';
% oripeplot.MarkerSize = markersizes;
% oripeplot.MarkerEdgeColor = colorset(3,:);
% 
% plot(templine,templine,'k:');
% legd = legend('mvSVDPE',"mvPE");
% legd.Box = "off";
% legd.Location = 'northwest';
% axis equal
% axis([50,100,50,100])
% hold off
% set(gca,'XTick',[50:10:100],'YTick',[50:10:100]);
% xlabel('10-fold cross-validation accuracy(%)','FontName','Times New Roman','FontSize',fonttxt_size,'FontWeight','bold')
% ylabel('Test accuracy(%)','FontName','Times New Roman','FontSize',fonttxt_size,'FontWeight','bold')
fig_description = '1_hf_mv_svdpe_pe';
lowerfont = 1;
box on

plotandprint;





% barplot = bar(y,'FaceColor','flat');
% colorset = hsv(6);
% for i = 1:4
%         barplot(i).FaceColor = colorset(i,:);
% 
%     % for j = 1:33
%     %barplot(i).CData(j,:) = colorset(i,:);
%     % end
% end