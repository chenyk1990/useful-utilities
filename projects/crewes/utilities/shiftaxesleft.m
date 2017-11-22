function shiftaxesleft(leftshift)

pos=get(gca,'position');

set(gca,'position',[pos(1)-leftshift pos(2:4)])