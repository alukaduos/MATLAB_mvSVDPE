clc;clear;close all;
load resulttable_080425.mat
num_classifer = size(resultsTablesvdpe,1);
x = resultsTablesvdpe.ModelType;
y = [resultsTablesvdpe.AccuracyValidation,resultsTablesvdpe.AccuracyTest,resultsTableoripe.AccuracyValidation, resultsTableoripe.AccuracyTest];
colorset = hsv(6); markersizes = 4;
templine = 50:0.1:100;
countf = 9;countdrawout = 10;
nrows = 2;ncol = 2;
posimat = figposi(nrows,ncol);linewidth_t = 1;
fonttxt_size = 8;fig_width = 17.6;fig_height = fig_width*0.5;
drawout = 1;
figure
hold on
svdpeplot = plot(y(:,1),y(:,2));
%svdpeplot.Color = colorset(1,:);
svdpeplot.Marker = 'o';
svdpeplot.LineStyle = 'none';
svdpeplot.MarkerSize = markersizes;
svdpeplot.MarkerEdgeColor = colorset(1,:);
oripeplot = plot(y(:,3),y(:,4));
oripeplot.Marker = '*';
oripeplot.LineStyle = 'none';
oripeplot.MarkerSize = markersizes;
oripeplot.MarkerEdgeColor = colorset(3,:);

plot(templine,templine,'k:');
legd = legend('H_{ES}',"H_P");
legd.Box = "off";
legd.Location = 'northwest';
axis equal
axis([50,100,50,100])
hold off
set(gca,'XTick',[50:10:100],'YTick',[50:10:100]);
xlabel('10-fold cross-validation accuracy(%)','FontName','Times New Roman','FontSize',fonttxt_size,'FontWeight','bold')
ylabel('Test accuracy(%)','FontName','Times New Roman','FontSize',fonttxt_size,'FontWeight','bold')
fig_description = '1_mv_svdpe_pe_gear';
lowerfont = 0;
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