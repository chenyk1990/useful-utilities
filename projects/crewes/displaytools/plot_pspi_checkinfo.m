function plot_pspi_checkinfo(checkinfo,thistitle)

exzos=checkinfo{1};
pmigs=checkinfo{2};
zchecks=checkinfo{3};

[nt,nx]=size(exzos{1});
[nz,nx2]=size(pmigs{1});

if(nx~=nx2)
    error(' x sizes of exzos and pmigs not the same');
end

nr=max([nt,nz]);

%determine overall scale factors
nsnaps=length(zchecks);
Amig=0;
Aex=0;
for k=1:nsnaps
    A=max(abs(pmigs{k}(:)));
    if(A>Amig);Amig=A;end
    A=max(abs(exzos{k}(:)));
    if(A>Aex);Aex=A;end
end

%make big pictures
bigpics=cell(size(zchecks));
titles=bigpics;
gap=10;
for k=1:length(bigpics)
    tmp=zeros(nr,2*nx+gap);
    tmp(1:nz,nx+gap+1:end)=pmigs{k}/Amig;
    tmp(1:nt,1:nx)=exzos{k}/Aex;
    bigpics{k}=tmp;
    titles{k}=['Extrapolated to depth ' int2str(zchecks(k))];
end

plotgathers(bigpics,nan,nan,'','',titles,thistitle,{{'Extrapolated zos''s',5,nr-10,10,'bold','left'},...
    {'Partial migrations',2*nx+gap,nr-10,10,'bold','right'}});

set(gcf,'name',['Post-stack PSPI viewer' thistitle]) 