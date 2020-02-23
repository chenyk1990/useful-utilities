%EMD based dip filter
clear;clc;close all;
% nx=256;dx=1;
% nt=512;dt=0.004;
% verb=0;
% 
% d0=zeros(512,256);
% d1=zeros(512,256);
% rsf_read(d0,'plane-c.rsf');
% rsf_read(d1,'plane.rsf');


%% make model 
% parameters definition
 flow=5;
 fhigh=120;
 seed=201314;
 nt=512;dt=0.004;
 nx=256;dx=1;
 f0=40;
 verb=0;
 
% making events
% %d1  =    linear_events(0.004,40,2,[0:10:10*79],1,0,1,2,2);
% d01 =   linear_events(dt,f0,(nt-1)*0.004,1:nx,1.2,0,1,200,2,'default');
% % figure;imagesc(d1);
% d02 =   linear_events(dt,f0,(nt-1)*0.004,1:nx,0.3,0.002,0.5,200,2,'default');
% % figure;imagesc(d2);
% d03 =   linear_events(dt,f0,(nt-1)*0.004,1:nx,0.5,0.005,0.5,200,2,'default');
% % figure;imagesc(d3);


%d1  =    linear_events(0.004,40,2,[0:10:10*79],1,0,1,2,2);
d01 =   linear_events(dt,f0,(nt-1)*0.004,1:nx,1.0,0,1,200,2,'default');
% figure;imagesc(d1);
d02 =   linear_events(dt,f0,(nt-1)*0.004,1:nx,0.6,0.002,0.5,200,2,'default');
% figure;imagesc(d2);
d03 =   linear_events(dt,f0,(nt-1)*0.004,1:nx,0.5,0.005,0.5,200,2,'default');
% figure;imagesc(d3);


d0=d01+d02+d03;d0=d0/max(max(d0));
randn('state',seed);
d1=d0+0.05*randn(nt,nx);
figure;imagesc([d0,d1]);

dip1=d1-fxemd(d1,5, 120, 0.004, 1, verb);
dip2=d1-fxemd(d1,5, 120, 0.004, 2, verb)-dip1;
dip3=d1-fxemd(d1,5, 120, 0.004, 3, verb)-dip1-dip2;
dip4=d1-fxemd(d1,5, 120, 0.004, 4, verb)-dip1-dip2-dip3;
dip5=d1-fxemd(d1,5, 120, 0.004, 5, verb)-dip1-dip2-dip3-dip4;

res=d1-dip1-dip2-dip3-dip4-dip5;


figure;imagesc([dip1,dip2,dip3;[dip4,dip5,res]]);


plane1=dip1+dip2;
plane2=dip3+dip4;
plane3=dip5+res;

figure;imagesc([plane1,plane2,plane3]);

%figure;imagesc([d1,dip1,dip2;[dip3,dip4,res]]);