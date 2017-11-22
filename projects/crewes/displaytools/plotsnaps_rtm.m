function plotsnaps_rtm(recsnaps,sousnaps,refsnaps,vel,xm,zm,times,modeltitle,units,cmap,cliplvls)
% PLOTSNAPS_RTM: works in conjunction with rtm_shot to display snapshots interactively
%
% plotsnaps_rtm(recsnaps,sousnaps,refsnaps,vel,xm,zm,times,modeltitle,units,cmap,cliplvls)
%
% recsnaps ... cell array of receiver snapshots from rtm_shot
% sousnaps ... cell array of source snapshots from rtm_shot
% refsnaps ... cell array of reflectivity snapshots from rtm_shot
% vel ... velocity model used in rtm_shot
% xm ... x coordinates for the velocity model
% zm ... z coordinates for the velocity model
% times ... vectors of times for each snapshot
% modeltitle ... string giving a title for the experiment
% units ... either 'metric' or 'ft' used for annotation
% ********* default 'metric' **********
% cmap ... colormap to use
% ********** default 'jet' ***********
% cliplvle ... vector of 3 clip levels for the initial display
% ******* default [3 3 8] ***********

if(nargin<11)
    cliplvls=[3 3 8];
end
if(nargin<10)
    cmap='jet';
end

if(nargin<9)
    units='metric';
end

if(strcmp(units,'metric'))
    xname='distance (m)';
    zname='depth (m)';
else
    xname='distance (ft)';
    zname='depth (ft)';
end


times=sort(times,'descend');
titles=cell(size(recsnaps));
for k=1:length(times)
    titles{k}=['time = ' time2str(times(k))];
end

ssize=get(0,'screensize');

plotsnaps(recsnaps,vel,xm,zm,xname,zname,titles,modeltitle,cmap,'Receiver data',6);
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
set(gcf,'position',[xnow,ynow,wid,ht],'name',['Receiver data snapshots ' modeltitle]);
bigfont(gcf,.75,1);
titlefontsize(.75,1);
plotsnaps(sousnaps,vel,xm,zm,xname,zname,titles,modeltitle,cmap,'Source model',6);
plotsnaps('group',hfig1);
hclip=findobj(gcf,'tag','clip');
hlevels=findobj(gcf,'tag','cliplevels');
levels=get(hlevels,'userdata');
ilevel=near(levels,cliplvls(2));
set(hclip,'value',ilevel(1));
plotsnaps('clip');
xnow=xnow+1.05*wid;
set(gcf,'position',[xnow,ynow,wid,ht],'name',['Source data snapshots ' modeltitle]);
bigfont(gcf,.75,1);
titlefontsize(.75,1);
plotsnaps(refsnaps,vel,xm,zm,xname,zname,titles,modeltitle,cmap,'Reflectivity',6);
plotsnaps('group',hfig1);
hclip=findobj(gcf,'tag','clip');
hlevels=findobj(gcf,'tag','cliplevels');
levels=get(hlevels,'userdata');
ilevel=near(levels,cliplvls(3));
set(hclip,'value',ilevel(1));
plotsnaps('clip');
ynow=ynow-1.2*ht;
xnow=xnow-.525*wid;
set(gcf,'position',[xnow,ynow,wid,ht],'name',['Reflectivity snapshots ' modeltitle]);
bigfont(gcf,.75,1);
titlefontsize(.75,1);