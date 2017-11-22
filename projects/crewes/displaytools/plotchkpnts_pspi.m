function plotchkpnts_pspi(chkpnts,t,vel,xm,zm,modeltitle,units,cmap,cliplvls)
% PLOTCHKPNTS_PSPI: works in conjunction with pspi_shot to display checkpoints interactively
%
% plotchkpnts_pspi(chkpnts,t,vel,xm,zm,modeltitle,units,cmap,cliplvls)
%
% chkpnts ... cell array of checkpoint info from pspi_shot
% vel ... velocity model used in pspi_shot
% xm ... x coordinates for the velocity model
% zm ... z coordinates for the velocity model
% modeltitle ... string giving a title for the experiment
% units ... either 'metric' or 'ft' used for annotation
% ********* default 'metric' **********
% cmap ... colormap to use
% ********** default 'jet' ***********
% cliplvle ... vector of 4 clip levels for the initial display
% ******* default [8 8 8 8] ***********

if(nargin<9)
    cliplvls=[8 8 8 8];
end
if(nargin<8)
    cmap='jet';
end

if(nargin<7)
    units='metric';
end

if(strcmp(units,'metric'))
    xname='distance (m)';
    zname='depth (m)';
else
    xname='distance (ft)';
    zname='depth (ft)';
end

depths=chkpnts{1};
reccheck=chkpnts{2};
soucheck=chkpnts{3};
cccheck=chkpnts{4};
deccheck=chkpnts{5};

depths=sort(depths);
titles=cell(size(reccheck));
for k=1:length(depths)
    titles{k}=['depth = ' num2str(depths(k))];
end

ssize=get(0,'screensize');
%bkgrnd=zeros(size(reccheck{1}));
tmp=linspace(0,1,length(t))';
bkgrnd=tmp(:,ones(size(xm)));
plotsnaps(reccheck,bkgrnd,xm,t,xname,zname,titles,modeltitle,cmap,'Receiver data',6);
hfig1=gcf;
plotsnaps('groupnew');
hclip=findobj(gcf,'tag','clip');
hlevels=findobj(gcf,'tag','cliplevels');
levels=get(hlevels,'userdata');
ilevel=near(levels,cliplvls(1));
set(hclip,'value',ilevel(1));
plotsnaps('clip');
xnow=ssize(3)/10;
ynow=ssize(4)/2;
ht=ssize(4)/3;
wid=ssize(3)/2.5;
set(gcf,'position',[xnow,ynow,wid,ht]);
bigfont(gcf,.75,1);
titlefontsize(.75,1);
plotsnaps(soucheck,bkgrnd,xm,t,xname,zname,titles,modeltitle,cmap,'Source model',6);
plotsnaps('group',hfig1);
hclip=findobj(gcf,'tag','clip');
hlevels=findobj(gcf,'tag','cliplevels');
levels=get(hlevels,'userdata');
ilevel=near(levels,cliplvls(2));
set(hclip,'value',ilevel(1));
plotsnaps('clip');
xnow=xnow+1.05*wid;
set(gcf,'position',[xnow,ynow,wid,ht]);
bigfont(gcf,.75,1);
titlefontsize(.75,1);
plotsnaps(cccheck,vel,xm,zm,xname,zname,titles,modeltitle,cmap,'Rcc',6);
plotsnaps('group',hfig1);
hclip=findobj(gcf,'tag','clip');
hlevels=findobj(gcf,'tag','cliplevels');
levels=get(hlevels,'userdata');
ilevel=near(levels,cliplvls(3));
set(hclip,'value',ilevel(1));
plotsnaps('clip');
ynow=ynow-1.2*ht;
xnow=xnow-1.05*wid;
set(gcf,'position',[xnow,ynow,wid,ht]);
bigfont(gcf,.75,1);
titlefontsize(.75,1);
plotsnaps(deccheck,vel,xm,zm,xname,zname,titles,modeltitle,cmap,'Rdec',6);
plotsnaps('group',hfig1);
hclip=findobj(gcf,'tag','clip');
hlevels=findobj(gcf,'tag','cliplevels');
levels=get(hlevels,'userdata');
ilevel=near(levels,cliplvls(4));
set(hclip,'value',ilevel(1));
plotsnaps('clip');
xnow=xnow+1.05*wid;
set(gcf,'position',[xnow,ynow,wid,ht]);
bigfont(gcf,.75,1);
titlefontsize(.75,1);