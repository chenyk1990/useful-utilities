function [ccx,xlags,ccy,ylags,xaqgrid,yaqgrid]=ccfoot(slice,x,y,dx,dy)
% CCFOOT ... measure footprint on a time slice
% [ccx,xlags,ccy,ylags,xgrid,ygrid]=ccfoot(slice,x,y,dx,dy)
%
% 
% slice ... input time slice
% x ... x (column) coordinate for slice
% y ... y (row) coordinate for slice
% dx ... line spacing in x during acquisition
% dy ... line spacing in y during acquisition
% ccx ... crosscorrelation in x between xgrid and slice
% xlags ... the xlags for ccx
% ccy ... crosscorrelation in y between ygrid and slice
% ylags ... the ylags for slice
% xgrid ... the xline acquisition grid
% ygrid the yline acquisition grid
%

%build xlinegrid
dxrow=abs(x(2)-x(1));
nxgrid=round(dx/dxrow);
xlinegrid=zeros(length(y),length(x)+2*nxgrid,'like',slice);
x2=(1:size(xlinegrid,2));

%build ylinegrid
dycol=abs(y(2)-y(1));
nygrid=round(dy/dycol);
ylinegrid=zeros(length(y)+2*nygrid,length(x),'like',slice);
y2=(1:size(ylinegrid,1))';

%populate grids
ix=nxgrid:nxgrid:length(x2)-1;
for k=1:length(ix)
    xlinegrid(:,ix(k))=1;
end
% seisplot(xlinegrid,y,x2);
iy=nygrid:nygrid:length(y2)-1;
for k=1:length(iy)
    ylinegrid(iy(k),:)=1;
end
% seisplot(ylinegrid,y2,x);

%bandlimit the grids
sigmax=.25;
sigmay=.25;
ynom=1:length(y);
xnom=1:length(x);
xlinegrid=wavenumber_gaussmask2(xlinegrid,sigmax,10);
xlinegrid=xlinegrid/max(xlinegrid(:));
ylinegrid=wavenumber_gaussmask2(ylinegrid,10,sigmay);
ylinegrid=ylinegrid/max(ylinegrid(:));
% figure
% subplot(2,1,1)
% seisplota(xlinegrid,ynom,x2);
% subplot(2,1,2)
% seisplota(ylinegrid,y2,xnom)
% prepfiga
% figure
% subplot(2,1,1)
% plot(x2,xlinegrid(round(length(ynom)/2),:))
% subplot(2,1,2)
% plot(y2,ylinegrid(:,round(length(xnom)/2)))

%need zero lag auto's for normalization
nx=length(x);
ny=length(y);
A0=sum(slice(:).^2);
Ax=sum(sum(xlinegrid(:,nxgrid+1:nxgrid+nx-1).^2));
Ay=sum(sum(ylinegrid(nygrid+1:nygrid+ny-1,:).^2));
%calculate the x correlation
xlags=-nxgrid:nxgrid;
nxlags=length(xlags);
ccx=zeros(1,nxlags);

% A=max(abs(slice(:)));
for k=1:nxlags
    x1=nxgrid+1+xlags(k);
    tmp=abs(slice).*xlinegrid(:,x1:x1+nx-1);
    %tmp=slice.*xlinegrid(:,x1:x1+nx-1);
    ccx(k)=sum(tmp(:))/sqrt(A0*Ax);
end
%calculate the y correlation
ylags=-nygrid:nygrid;
nylags=length(ylags);
ccy=zeros(1,nylags);

for k=1:nylags
    y1=nygrid+1+ylags(k);
    tmp=abs(slice).*ylinegrid(y1:y1+ny-1,:);
    %tmp=slice.*ylinegrid(y1:y1+ny-1,:);
    ccy(k)=sum(tmp(:))/sqrt(A0*Ay);
end
% figure
% subplot(1,2,1)
% plot(xlags,ccx)
% subplot(1,2,2)
% plot(ylags,ccy)

xaqgrid=xlinegrid(:,nxgrid+1:nxgrid+nx-1);
yaqgrid=ylinegrid(nygrid+1:nygrid+ny-1,:);

