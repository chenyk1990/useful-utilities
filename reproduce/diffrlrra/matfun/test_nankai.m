clc;close all;clear;

% d=zeros(1750,401);
% rsf_read(d,'../nankai/vc15.rsf');
% 
% n1=zeros(1750,401);
% rsf_read(n1,'../nankai/pwd-vc.rsf');

% save nankai.mat d n1
addpath(genpath('/Users/chenyk/chenyk/matlibcyk'));
load nankai.mat

n1=-n1;
d1=d-n1;

figure;imagesc([d,d1,n1]);colormap(seis);caxis([-100,100]);


d2=fxmssa_win(d,0,120,0.004,5,0,200,50,0.5,0.5);
n2=d-d2;
figure;imagesc([d,d2,n2]);colormap(seis);caxis([-100,100]);

d3=fxdmssa_win(d,0,120,0.004,5,3,0,200,50,0.5,0.5);
n3=d-d3;
figure;imagesc([d,d3,n3]);colormap(seis);caxis([-100,100]);
%d3 and d33 exactly the same
d33=fxydmssa_win(d,0,120,0.004,2,3,0,100,20,1,0.5,0.5,0.5);
n33=d-d33;
figure;imagesc([d,d33,n33]);colormap(seis);caxis([-100,100]);

% automatic rank
d333=fxydmssa_win_auto(d,0,120,0.004,6,3,0,100,20,1,0.5,0.5,0.5,2);
n333=d-d333;
figure;imagesc([d,d333,n333]);colormap(seis);caxis([-100,100]);

d222=fxymssa_win_auto(d,0,120,0.004,6,0,100,20,1,0.5,0.5,0.5,2);
n222=d-d222;
figure;imagesc([d,d222,n222]);colormap(seis);caxis([-100,100]);


%%
figure;imagesc([n1,n2,n3]);colormap(seis);caxis([-100,100]);


test1=d-fxymssa_win(d,0,120,0.004,3,0,100,20,1,0.5,0.5,0.5);
test2=d-fxydmssa_win(d,0,120,0.004,3,3,0,100,20,1,0.5,0.5,0.5);
figure;imagesc([test1,test2]);colormap(seis);caxis([-100,100]);



