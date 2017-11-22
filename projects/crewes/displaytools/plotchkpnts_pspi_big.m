function plotchkpnts_pspi_big(chkpnts,vel,modeltitle,cmap,gfacts,alevel)


if(nargin<6)
    alevel=1;
end
if(nargin<5)
    gfacts=[1 1 1 1];
end

if(nargin<4)
    cmap='jet';
end

depths=chkpnts{1};
reccheck=chkpnts{2};
soucheck=chkpnts{3};
cccheck=chkpnts{4};
deccheck=chkpnts{5};

%determine size
[nt,nx1]=size(reccheck{1});
[nz,nx2]=size(cccheck{1});
ny=max([nt nz]);
nx=max([nx1 nx2]);
ngx=round(.1*nx);%gap
ngy=round(.1*ny);%gap
nxbig=2*nx+ngx;
nybig=2*ny+ngx;

%determine overall scale factors
nsnaps=length(reccheck);
Arec=0;
Asou=0;
Acc=0;
Adec=0;
for k=1:nsnaps
    A=max(abs(reccheck{k}(:)));
    if(A>Arec);Arec=A;end
    A=max(abs(soucheck{k}(:)));
    if(A>Asou);Asou=A;end
    A=max(abs(cccheck{k}(:)));
    if(A>Acc);Acc=A;end
    A=max(abs(deccheck{k}(:)));
    if(A>Adec);Adec=A;end
end

%build big snapshots
bigsnaps=cell(1,nsnaps);

for k=1:nsnaps
    tmp=zeros(nybig,nxbig);
    tmp(1:nt,1:nx1)=reccheck{k}*gfacts(1)/Arec;
    tmp(1:nt,nx1+ngx+1:2*nx1+ngx)=soucheck{k}*gfacts(2)/Asou;
    tmp(nt+ngy+1:nt+ngy+nz,1:nx2)=cccheck{k}*gfacts(3)/Acc;
    tmp(nt+ngy+1:nt+ngy+nz,nx2+ngx+1:nx2+ngx+nx2)=deccheck{k}*gfacts(4)/Adec;
    bigsnaps{k}=tmp;
end
%build bigvel
tmp=linspace(min(vel(:)),max(vel(:)),nt)';
bigvel=zeros(nybig,nxbig)+min(vel(:));
bigvel(1:nt,1:nx1)=tmp(:,ones(1,nx1));
bigvel(1:nt,nx1+ngx+1:2*nx2+ngx)=tmp(:,ones(1,nx1));
bigvel(nt+ngy+1:nt+ngy+nz,1:nx2)=vel;
bigvel(nt+ngy+1:nt+ngy+nz,nx2+ngx+1:nx2+ngx+nx2)=vel;

%depths=sort(depths);
titles=cell(size(reccheck));
for k=1:length(depths)
    titles{k}=['depth = ' num2str(depths(k))];
end

plotsnaps(bigsnaps,bigvel,nan,nan,'','',titles,modeltitle,cmap,'PSPI shot record migration',10,...
    {{'Receiver data',round(ngx/5),round(ngy/3),10,'bold','left'},{'Source model',nxbig-round(ngx/5),round(ngy/3),10,'bold','right'}...
    ,{'Ref CCIM',round(ngx/5),nybig-round(ngy/2),10,'bold','left'},{'Ref DECIC',nxbig-round(ngx/5),nybig-round(ngy/2),10,'bold','right'}});
hclip=findobj('tag','clip');
hlevels=findobj(gcf,'tag','cliplevels');
levels=get(hlevels,'userdata');
ilevel=near(abs(levels),alevel);
set(hclip,'value',ilevel(1));
plotsnaps('clip');
set(gcf,'name','PLOTSNAPS_BIG: PSPI shot record migration')