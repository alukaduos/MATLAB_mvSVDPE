savefigpath = 'D:\MATLAB\Work\MTS_Entropy\figureset\';

set(gcf,'Position',posimat(mod(countf,nrows*ncol)+1,:));
set(gcf,'PaperUnits','centimeters','PaperPosition',[5,5,fig_width,fig_height]);
countf = countf+1;
if lowerfont
    set(gca,"FontSize",fonttxt_size-3)
else
    set(gca,"FontSize",fonttxt_size)
end
  if lowerfont
    tempt1.FontSize = fonttxt_size;
  end
try
set(gca,"FontName",'Times New Roman','FontWeight','bold')
%set(gca,"FontName",'Times New Roman','FontWeight','bold')

catch
% set(gca,"FontSize",fonttxt_size)
   
end
picname = strcat(savefigpath,'fig_',num2str(countdrawout),'_',fig_description);

if drawout
    print(picname,'-depsc','-r300','-vector')
    print(picname,'-dpng','-r300','-vector')
    countdrawout = countdrawout+1;
    drawout = 0;
end
