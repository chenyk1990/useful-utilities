function haxes=tvphasepro(seis,t,x,xref,t1s,xs,twins,pos,klrs,delx)
% TVPHASEPRO: compute and display relative phase profiles in 3 time windows
%
% haxes=tvphasepro(seis,t,x,xref,t1s,xs,twins,pos,klrs,delx)
%
% seis ... seismic section
% t ... time coordinate for seis
% x ... space coordinate for seis
% xref ... reference trace location in x units
% t1s ... matrix of 3 rows by n columns (n may be 1) giving the window center times at the x
%       coordinates in xs for 3 windows. Number of columns is length(xs). Number of rows is always
%       3.
% xs ... vector of x coordinates at which the t1s are specified. length(xs) must equal size(t1s,2)
% twins ... vector of length 3 giving the temporal window sizes for the three phase analysis windows
%       May be a scalar if all windows are the same
% pos ... length 4 vector of [left, bottom, width, height] specifies the location and size of the 
%     box that the three xes will appear in, relative to the lower-left corner of the Figure window,  
%     in normalized units where (0,0) is the lower-left corner and (1.0,1.0) is the upper-right.
%     Can be specified as an axes handle. In which case, pos is take from the handle, then the axes
%     is deleted and replaed by 3 samaller axes.
% klrs ... 3x3 matrix of colors (rgb) to plot with. Window 1 is klrs(1,:) etc.
% ************ default uses the standard color order of the axes ***************
% delx ... half width of the trace ensemble size that is compared to the reference trace. See
%     phaseprofile.m for more information. This is specficied as a number of traces not in x units.
% ************ default = 20 (traces) ***********
%
% haxes ... [haxp haxd haxcc] three axes handles for phase delay and cc

nxs=length(xs);
[mt,nt]=size(t1s);
if(mt~=3)
    error('t1s must have 3 rows, one per window');
end
if(nt~=nxs)
    error('length(xs) must equal size(t1s,2)')
end
if(length(twins)==1)
    twins=twins*ones(1,3);
end
if(nargin<9)
    klrs=get(gca,'colororder');
end
if(nargin<10)
    delx=20;
end

[nt,nx]=size(seis);
if(length(t)~=nt)
    error('t and seis have incompatible sizes');
end
if(length(x)~=nx)
    error('x and seis have incompatible sizes');
end
if(isgraphics(pos))
    haxes=pos;
    pos=get(haxes,'position');
    delete(haxes);
end

phs=zeros(3,nx);
delay=phs;
cc=phs;

for k=1:3
    [phs(k,:),delay(k,:),cc(k,:)]=phaseprofile(seis,t,x,xref,t1s(k,:),xs,twins(k),delx);
end

xmin=min(x);
xmax=max(x);
xm=.5*(xmin+xmax);

ht=pos(4)/3;
haxp=axes('position',[pos(1) pos(2)+2*ht pos(3) ht]);
hh=plot(x,phs);
set(hh(1),'color',klrs(1,:));set(hh(2),'color',klrs(2,:));set(hh(3),'color',klrs(3,:));
set(haxp,'xticklabel','','xgrid','on','ygrid','on');
ylabel('degrees');
xlim([xmin xmax]);
ylim([-200 200])
yl=get(gca,'ylim');
dy=abs(diff(yl));
line([xref xref],yl,'color','k','linestyle','--');
ixref=near(x,xref);
patch([x(ixref-delx) x(ixref+delx) x(ixref+delx) x(ixref-delx)],[yl(2) yl(2) yl(1) yl(1)],-1*ones(1,4),...
    .9*ones(1,3),'edgecolor',.9*ones(1,3),'tag','ensemble');
text(xref,yl(2)+.1*dy,'Reference location','horizontalalignment','center')
title('phase','position',[xm yl(1)+.25*dy -eps])

haxd=axes('position',[pos(1) pos(2)+ht pos(3) ht]);
hh=plot(x,delay);
set(hh(1),'color',klrs(1,:));set(hh(2),'color',klrs(2,:));set(hh(3),'color',klrs(3,:));
set(haxd,'xticklabel','','xgrid','on','ygrid','on','yaxislocation','right');
ylabel('seconds')
xlim([xmin xmax]);
yl=get(gca,'ylim');
dy=abs(diff(yl));
line([xref xref],yl,'color','k','linestyle','--');
title('delay','position',[xm yl(1)+.25*dy -eps])

haxcc=axes('position',[pos(1:3) ht]);
hh=plot(x,cc);
set(hh(1),'color',klrs(1,:));set(hh(2),'color',klrs(2,:));set(hh(3),'color',klrs(3,:));
xlim([xmin xmax]);
yl=get(gca,'ylim');
dy=abs(diff(yl));
title('CC','position',[xm yl(1)+.25*dy -eps])
set(haxcc,'xgrid','on','ygrid','on');
ylabel('cc value')
line([xref xref],yl,'color','k','linestyle','--');
xlim([min(x) max(x)]);

haxes=[haxp haxd haxcc];
