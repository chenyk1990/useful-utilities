function [seisout,opw,xop,top]=depthstep(seis,t,x,v,dz,w,tw,ieveryt,ieveryx,dipmax)

if(nargin<10)
    dipmax=70;
end
if(nargin<9)
    ieveryx=nan;
end
if(nargin<8)
    ieveryt=nan;
end

if(isnan(ieveryx))
    ieveryx=1;
end
if(isnan(ieveryt))
    ieveryt=1;
end
nt=length(t);
nx=length(x);
pctx=10;
pctt=10;
%make the operator
tnot=2*dz/v;
xnot=x(round(nx/2));
dips=atan2d((x-xnot),dz);
indx=find(abs(dips)<=dipmax);
tmax=1.2*2*dz/(v*cosd(dipmax));
indt=near(t,0,tmax)';
op=zeros(length(indt),length(indx));
xop=x(indx)-xnot;
top=t(indt);
op=event_hyp(op,top,xop,tnot,0,v,1,3);
%apply tapers with mwindow
mwx=mwindow(length(indx),pctx)';
mwt=mwindow(length(indt),pctt);
op=op.*mwx(ones(size(indt)),:).*mwt(:,ones(size(indx)));
%apply wavelet
nzero=near(tw,0);
opw=flipud(convz(op,w,nzero));


%now the convolution

indt=1:ieveryt:nt;
indx=1:ieveryx:nx;
tmp=zeros(size(seis));
tmp(indt,indx)=seis(indt,indx);
tmp2=conv2(tmp,opw,'same');

seisout=tmp2*norm(tmp(:))/norm(tmp2(:));


