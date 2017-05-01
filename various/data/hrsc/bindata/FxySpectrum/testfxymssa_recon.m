%  Simultaneous seismic data denoising and reconstruction via multichannel singular spectrum analysis

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
%  Reference:   Simultaneous seismic data denoising and reconstruction via multichannel
%               singular spectrum analysis, Geophysics, 2011, 76, V25-V32
%

clc;clear;close all;
%% generate synthetic data
nt=300;
nx=20;
ny=20;
nf=2^nextpow2(nt);


fp=fopen('origin.bin','rb');
d=fread(fp,[nt inf],'float');
d=reshape(d,nt,nx,ny);
fclose(fp);

fp=fopen('data_with_noise.bin','rb');
dn=fread(fp,[nt inf],'float');
dn=reshape(dn,nt,nx,ny);
fclose(fp);

fp=fopen('decimate.bin','rb');
decimate=fread(fp,[nt inf],'float');
decimate=reshape(decimate,nt,nx,ny);
fclose(fp);

fp=fopen('mssa.bin','rb');
mssa=fread(fp,[nt inf],'float');
mssa=reshape(mssa,nt,nx,ny);
fclose(fp);

fp=fopen('rs.bin','rb');
rs=fread(fp,[nt inf],'float');
rs=reshape(rs,nt,nx,ny);
fclose(fp);

% Transform into F-X domain
DATA_FX_d=fft(d,nf,1);
DATA_FX_dn=fft(dn,nf,1);
DATA_FX_decimate=fft(decimate,nf,1);
DATA_FX_mssa=fft(mssa,nf,1);
DATA_FX_rs=fft(rs,nf,1);

% Extracing FXY Spectrum
nff=11;
Spec_d=abs(squeeze(DATA_FX_d(nff,:,:)));
Spec_dn=abs(squeeze(DATA_FX_dn(nff,:,:)));
Spec_decimate=abs(squeeze(DATA_FX_decimate(nff,:,:)));
Spec_mssa=abs(squeeze(DATA_FX_mssa(nff,:,:)));
Spec_rs=abs(squeeze(DATA_FX_rs(nff,:,:)));

fp=fopen('FXYspec_clean.bin','wb');
fwrite(fp,Spec_d,'float');
fclose(fp);

fp=fopen('FXYspec_noisy.bin','wb');
fwrite(fp,Spec_dn,'float');
fclose(fp);

fp=fopen('FXYspec_decimate.bin','wb');
fwrite(fp,Spec_decimate,'float');
fclose(fp);

fp=fopen('FXYspec_mssa.bin','wb');
fwrite(fp,Spec_mssa,'float');
fclose(fp);

fp=fopen('FXYspec_hrsc.bin','wb');
fwrite(fp,Spec_rs,'float');
fclose(fp);



% Spec_d=abs(squeeze(DATA_FX_d(:,nff,:)));
% Spec_dn=abs(squeeze(DATA_FX_dn(:,nff,:)));
% Spec_decimate=abs(squeeze(DATA_FX_decimate(:,nff,:)));
% Spec_mssa=abs(squeeze(DATA_FX_mssa(:,nff,:)));
% Spec_rs=abs(squeeze(DATA_FX_rs(:,nff,:)));


% Ploting
% colormap(jet);
% subplot(2,3,1);
% subimage(Spec_d);
% subplot(2,3,2);
% subimage(Spec_dn);
% subplot(2,3,3);
% subimage(Spec_decimate);
% subplot(2,3,4);
% subimage(Spec_mssa);
% subplot(2,3,5);
% subimage(Spec_rs);

% Ploting
figure(1)
colormap(jet);
subplot(2,3,1);
imagesc(Spec_d);
subplot(2,3,2);
imagesc(Spec_dn);
subplot(2,3,3);
imagesc(Spec_decimate);
subplot(2,3,4);
imagesc(Spec_mssa);
subplot(2,3,5);
imagesc(Spec_rs);

nff=10;
Spec_d=abs(squeeze(DATA_FX_d(1:120,nff,:)));
Spec_dn=abs(squeeze(DATA_FX_dn(1:120,nff,:)));
Spec_decimate=abs(squeeze(DATA_FX_decimate(1:120,nff,:)));
Spec_mssa=abs(squeeze(DATA_FX_mssa(1:120,nff,:)));
Spec_rs=abs(squeeze(DATA_FX_rs(1:120,nff,:)));

fp=fopen('FXspec_clean.bin','wb');
fwrite(fp,Spec_d,'float');
fclose(fp);

fp=fopen('FXspec_noisy.bin','wb');
fwrite(fp,Spec_dn,'float');
fclose(fp);

fp=fopen('FXspec_decimate.bin','wb');
fwrite(fp,Spec_decimate,'float');
fclose(fp);

fp=fopen('FXspec_mssa.bin','wb');
fwrite(fp,Spec_mssa,'float');
fclose(fp);

fp=fopen('FXspec_hrsc.bin','wb');
fwrite(fp,Spec_rs,'float');
fclose(fp);


figure(2)
colormap(jet);
subplot(2,3,1);
imagesc(Spec_d);
subplot(2,3,2);
imagesc(Spec_dn);
subplot(2,3,3);
imagesc(Spec_decimate);
subplot(2,3,4);
imagesc(Spec_mssa);
subplot(2,3,5);
imagesc(Spec_rs);

% figure(1);
% imagesc(Spec_d);
% figure(2);
% imagesc(Spec_dn);
% figure(3);
% imagesc(Spec_decimate);
% figure(4);
% imagesc(Spec_mssa);
% figure(5);
% imagesc(Spec_rs);


