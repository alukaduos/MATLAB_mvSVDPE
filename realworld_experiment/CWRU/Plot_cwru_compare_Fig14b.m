clc;clear;close all;
load result_cwru_090425.mat
%load test4_cwru_080225.mat
num_classifer = size(resultcwru25032070train,1);
x = resultcwru25032070train.Preset;
y = [resultcwru25032070train.ra1,resultcwru25032070train.rt1,resultcwru25032070train.ra2,resultcwru25032070train.rt2];
colorset = hsv(6); markersizes = 4;
templine = 70:0.1:100;
countf = 15;countdrawout = 15;
nrows = 2;ncol = 2;
posimat = figposi(nrows,ncol);linewidth_t = 1;
fonttxt_size = 8;fig_width = 8;fig_height = fig_width;
drawout = 1;
figure
hold on
svdpeplot = plot(y(:,1),y(:,2));
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
axis([70,105,70,105])
hold off
set(gca,'XTick',[70:5:100],'YTick',[70:5:100]);
xlabel('10-fold cross-validation accuracy(%)','FontName','Times New Roman','FontSize',fonttxt_size,'FontWeight','bold')
ylabel('Test accuracy(%)','FontName','Times New Roman','FontSize',fonttxt_size,'FontWeight','bold')
fig_description = '3_mv_svdpe_pe_cwru_enlarge';
lowerfont = 0;
box on
% countf = countf-1;countdrawout = countdrawout-1;
drawout = 1;
plotandprint;





