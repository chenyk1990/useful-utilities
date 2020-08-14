clc;clear;close all;

d=zeros(800,501);
rsf_read(d,'../synth/data.rsf');
figure;imagesc(d);

%% GRR
d1=fxymssa_win(d,0,120,0.004,5,0,800,501,1,0.5,0.5,0.5);
d11=fxymssa_win_auto(d,0,120,0.004,10,0,800,501,1,0.5,0.5,0.5,2);

figure;imagesc([d,d1,d-d1]);
figure;imagesc([d,d11,d-d11]);

%% LRR
d2=fxymssa_win(d,0,120,0.004,3,0,200,100,1,0.5,0.5,0.5);figure;imagesc([d,d2,d-d2]);
d22=fxymssa_win_auto(d,0,120,0.004,5,0,200,100,1,0.5,0.5,0.5,2);

figure;imagesc([d,d2,d-d2]);caxis([-0.5,0.5]);colormap(seis);
figure;imagesc([d,d22,d-d22]);caxis([-2,2]);colormap(seis);

%% PWD
n3=zeros(800,501);
rsf_read(n3,'../synth/diffr-pwd-n.rsf');
d3=d+n3;
figure;imagesc([d,d3,d-d3]);caxis([-0.5,0.5]);colormap(seis);

dip=zeros(800,501);
rsf_read(dip,'../synth/data-dip.rsf');
figure;imagesc([dip]);colormap(jet);

%