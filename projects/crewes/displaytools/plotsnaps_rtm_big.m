function plotsnaps_rtm_big(recsnaps,sousnaps,refsnaps,vel,xm,zm,times,modeltitle,units,cmap,gfacts)



if(nargin<11)
    gfacts=[1 1 1];
end

if(nargin<10)
    cmap='jet';
end

% if(nargin<9)
%     units='metric';
% end

% if(strcmp(units,'metric'))
%     xname='distance (m)';
%     zname='depth (m)';
% else
%     xname='distance (ft)';
%     zname='depth (ft)';
% end


times=sort(times,'descend');
titles=cell(size(recsnaps));
for k=1:length(times)
    titles{k}=['time = ' time2str(times(k))];
end

%determine size
[ny,nx]=size(recsnaps{1});
ngx=round(.1*nx);
ngy=round(.1*ny);
nx2=2*nx+ngx;
ny2=2*ny+ngx;

%determine overall scale factors
nsnaps=length(recsnaps);
Arec=0;
Asou=0;
Aref=0;
for k=1:nsnaps
    A=max(abs(recsnaps{k}(:)));
    if(A>Arec);Arec=A;end
    A=max(abs(sousnaps{k}(:)));
    if(A>Asou);Asou=A;end
    A=max(abs(refsnaps{k}(:)));
    if(A>Aref);Aref=A;end
end

%build big snapshots
bigsnaps=cell(1,nsnaps);
for k=1:nsnaps
    tmp=zeros(ny2,nx2);
    tmp(1:ny,1:nx)=recsnaps{k}*gfacts(1)/Arec;
    tmp(1:ny,nx+ngx+1:2*nx+ngx)=sousnaps{k}*gfacts(2)/Asou;
    tmp(ny+ngy+1:2*ny+ngy,round((nx+ngx)/2)+1:round((nx+ngx)/2)+nx)=refsnaps{k}*gfacts(3)/Aref;
    bigsnaps{k}=tmp;
end
bigvel=zeros(ny2,nx2);
bigvel(1:ny,1:nx)=vel;
bigvel(1:ny,nx+ngx+1:2*nx+ngx)=vel;
bigvel(ny+ngy+1:2*ny+ngy,round((nx+ngx)/2)+1:round((nx+ngx)/2)+nx)=vel;

plotsnaps(bigsnaps,bigvel,nan,nan,'','',titles,modeltitle,cmap,'RTM shot record migration',10,...
    {{'Receiver data',round(ngx/5),ny+round(ngy/3),10,'bold','left'},{'Source model',nx2-round(ngx/5),ny+round(ngy/3),10,'bold','right'}...
    ,{'Reflectivity',round(nx2/2),ny2-round(ngy/2),10,'bold','center'}});
hclip=findobj('tag','clip');
hlevels=findobj(gcf,'tag','cliplevels');
levels=get(hlevels,'userdata');
ilevel=near(levels,8);
set(hclip,'value',ilevel);
plotsnaps('clip');
set(gcf,'name','PLOTSNAPS_BIG: RTM shot record migration')