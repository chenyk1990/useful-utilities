function eqzoom

haxes=findobj(gcf,'type','axes');
xl=get(gca,'xlim');yl=get(gca,'ylim');
set(haxes,'xlim',xl,'ylim',yl);