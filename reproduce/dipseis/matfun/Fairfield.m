%fairfieldnodal

clear;clc;close all;
dt=0.004;
verb=1;
d1=zeros(3001,560);
rsf_read(d1,'uorig.rsf')
d1=d1(1:2000,1:560);
figure;imagesc(d1,[-0.001,0.001]);

%d2=fxdecon(d1,5,120,0.004,10,0.01,1);
%figure;imagesc([d1,d2,d1-d2],[-0.001,0.001]);

%d3=fxemd(d1,5,120,0.004,1,1);
%figure;imagesc([d1,d3,d1-d3],[-0.002,0.002]);

dip1=d1-fxemd(d1,5, 120, 0.004, 1, verb);
dip2=d1-fxemd(d1,5, 120, 0.004, 2, verb)-dip1;
dip3=d1-fxemd(d1,5, 120, 0.004, 3, verb)-dip1-dip2;
dip4=d1-fxemd(d1,5, 120, 0.004, 4, verb)-dip1-dip2-dip3;
res=d1-dip1-dip2-dip3-dip4;
figure;imagesc([[d1,dip1,dip2];[dip3,dip4,res]],[-0.002,0.002]);

dips=[d1,dip1,dip2,dip3,dip4,res];

rsf_create('dips.rsf',size(dips)')
rsf_write(dips,'dips.rsf');
