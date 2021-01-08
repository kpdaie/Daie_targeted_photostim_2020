% function figure_initialize
clf
f=gcf;
whitebg('w')
set(f,'color','w');
set(f,'units','inches');
% set(f,'defaultaxescolororder',[0 0 0])
% orient portrait
set(f,'paperposition',[1.25 2.25 6.5 9])
set(f,'renderer','painters')

% % set(f,'defaultaxeslinewidth',1)
% set(f,'defaultaxesfontweight','normal')
% set(f,'defaultaxesfontname','Arial')
% set(f,'defaultaxesfontsize',8)
% % set(f,'defaultaxesticklength',[.02 .025])
% % set(f,'defaultaxestickdir','out')
% set(f,'defaultaxescolor','none')
set(f,'defaultaxesunits','inches')
% cllL=get(f,'DefaultAxesColorOrder');
% cllL(2,:)=[1 0 0];cllL(3,:)=[.5 0 0];
% set(f,'DefaultAxesColorOrder',cllL);
% 
% set(f,'defaulttextfontweight','normal')
% set(f,'defaulttextfontname','Arial')
% set(f,'defaulttextfontsize',8)

% set(f,'defaultlinelinewidth',3)
set(f,'defaultlinecolor',[0 0 0])
set(gcf,'units','inches');


set(gcf, 'PaperPositionMode', 'auto');
