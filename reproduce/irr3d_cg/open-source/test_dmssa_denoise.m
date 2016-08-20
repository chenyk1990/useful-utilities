% Demonstration script for 
% seismic denoising via
% damped multichannel singular spectrum analysis
%
%  Copyright (C) 2015 The University of Texas at Austin
%  Copyright (C) 2015 Yangkang Chen
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%
%  Reference:   
%  [1] Oropeza, V., and M. Sacchi, 2011, Simultaneous seismic data denoising and reconstruction via multichannel singular spectrum analysis, Geophysics, 76, V25-V32.
%  [2] Huang, W., R. Wang, M. Zhang, and Y. Chen, 2015, Damped multichannel singular spectrum analysis for 3D random noise attenuation: SEG expanded abstracts: 85th Annual international meeting, 4714–4719.
%  [3] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, XX, XX-XX.
%				
%

clc;clear;close all;

%% generate synthetic data
a1=zeros(300,20);
[n,m]=size(a1);
a3=a1;
a4=a1;

k=0;
a=0.1;
b=1;
for t=-0.055:0.002:0.055
    k=k+1;
    b1(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
    b2(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b3(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b4(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
end
for i=1:m
  t1(i)=round(140);
  t3(i)=round(-6*i+180);
  t4(i)=round(6*i+10);
  a1(t1(i):t1(i)+k-1,i)=b1; 
  a3(t3(i):t3(i)+k-1,i)=b1; 
  a4(t4(i):t4(i)+k-1,i)=b1;
end

temp=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
for j=1:20
    a4=zeros(300,20);
    for i=1:m
  t4(i)=round(6*i+10+3*j); 
  a4(t4(i):t4(i)+k-1,i)=b1;
  
  t1(i)=round(140-2*j);
  a1(t1(i):t1(i)+k-1,i)=b1;
    end
    shot(:,:,j)=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
end
plane3d=shot;
d=plane3d/max(max(max(plane3d)));

%% adding noise
randn('state',201314);
var=0.2;
dn=d+var*randn(size(d));

%% denoise
flow=0;fhigh=250;dt=0.004;N=3;verb=1;
d1=fxymssa(dn(:,:,:),flow,fhigh,dt,N,verb);
figure;imagesc([d(:,:,9),dn(:,:,9),d1(:,:,9),dn(:,:,9)-d1(:,:,9)]);caxis([-0.5,0.5]);

%% denoise
flow=0;fhigh=250;dt=0.004;N=3;verb=1;K=3
d2=fxydmssa(dn(:,:,:),flow,fhigh,dt,N,K,verb);
figure;imagesc([d(:,:,9),dn(:,:,9),d2(:,:,9),dn(:,:,9)-d2(:,:,9)]);caxis([-0.5,0.5]);



