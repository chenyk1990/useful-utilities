
clc;clear;close all;
%first scons in /Users/chenyk/chenyk/diffr/gom/
gom=zeros(301,920);
d0=zeros(301,920);
% rsf_read(gom,'/Users/chenyk/chenyk/diffr/gom/gom.rsf');
fid=fopen('/Users/chenyk/chenyk/diffr/gom/gom.bin','r');
gom=fread(fid,[301,920],'float');
d=gom;

fid=fopen('/Users/chenyk/chenyk/diffr/gom/gom-shp.bin','r');
d0=fread(fid,[301,920],'float');
% rsf_read(d0,'/Users/chenyk/chenyk/diffr/gom/gom-shp.rsf');

figure;imagesc([gom,gom-d0,d0]);colormap(seis);caxis([-2000,2000]);
figure;imagesc(d0);colormap(seis);caxis([-2000,2000]);

%% 
% d3=fxmssa(d,0,100,0.004,45,0);
d3=fxymssa(d,0,100,0.004,50,0);
figure;imagesc([d,d3,d-d3]);caxis([-2000,2000]);colormap(seis);


% not good actually (since no random noise)
d4=fxydmssa(d,0,100,0.004,45,5,0);
figure;imagesc([d,d4,d-d4]);caxis([-2000,2000]);colormap(seis);
figure;imagesc([d-d4,d0]);colormap(seis);caxis([-3000,3000]);


figure;imagesc([d-d3,d0]);colormap(seis);caxis([-3000,3000]);
annotation(gcf,'arrow',[0.469444444444445 0.457638888888889],...
    [0.608527131782946 0.566815953104883],'Color',[1 0 0],'LineWidth',4);

annotation(gcf,'arrow',[0.855555555555556 0.84375],...
    [0.606589147286821 0.564877968608759],'Color',[1 0 0],'LineWidth',4);

annotation(gcf,'arrow',[0.609722222222223 0.597916666666667],...
    [0.674418604651163 0.6327074259731],'Color',[1 0 0],'LineWidth',4);

annotation(gcf,'arrow',[0.218055555555556 0.206250000000001],...
    [0.674418604651162 0.6327074259731],'Color',[1 0 0],'LineWidth',4);

annotation(gcf,'arrow',[0.634722222222223 0.622916666666667],...
    [0.825581395348837 0.783870216670775],'Color',[1 0 0],'LineWidth',4);

annotation(gcf,'arrow',[0.254166666666668 0.242361111111112],...
    [0.823643410852712 0.78193223217465],'Color',[1 0 0],'LineWidth',4);


fid=fopen('/Users/chenyk/chenyk/diffr/gom/gom-rr.bin','w');
fwrite(fid,d3,'float');
