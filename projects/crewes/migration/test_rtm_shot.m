%% model a shot record over the wedge model
dx=5;
xmax=2500;zmax=1000; %maximum line length and maximum depth
vhigh=4000;vlow=2000; % high and low velocities
[vel,x,z]=wedgemodel(dx,xmax,zmax,vhigh,vlow);
fdom=35;
xshot=max(x)/3;
zshot=0.0;
dt=.004; %temporal sample rate
dtstep=.0005;%time step size
%[w,tw]=ricker(dtstep,fdom,.2);%wavelet
[w,tw]=wavemin(dtstep,fdom,.2);
tmax=2*zmax/vlow; %maximum time
xrec=x;
zrec=zeros(size(xrec));
[shot,t]=afd_shotrec_alt(dx,dtstep,dt,tmax,vel,xshot,zshot,xrec,zrec,w,tw,2);

seisplotfk(shot,t,x,'Test shot')
%% rtm_shot migration of the shot record from the previous section with exact velocity model
dtstep=.0005;
itcorr=4;
tout=0:.05:tmax;
xv=x;
lap=2;
illume=1;
[shotm,shotmf,xm,zm,recsnaps,sousnaps,refsnaps,xms,zms]=rtm_shot(shot,t,xrec,zrec,xshot,zshot,w,tw,vel,xv,...
    dtstep,itcorr,tout,lap,illume);

seisplottwo(shot,t,xrec,'Input shot record',shotm,zm,xm,['RTM migrated shot record (illume=' int2str(illume) ')'])

seisplottwo(shotm,zm,xm,'RTM migrated shot record',shotmf,zm,xm,'RTM migration with Laplacian filter')

plotsnaps_rtm(recsnaps,sousnaps,refsnaps,vel,xm,zm,tout,'wedgemodel','metric','jet')

%% rtm_shot of the shot record from the first section with smooth velocity model
%smooth the model by convolution with a Gaussian whose halfwidth is the dominant wavelength
vmean=mean(vel(:));
lambda=vmean/fdom;
velsmo=gaussian_smoother(vel,x,z,.1*lambda);
figure
subplot(1,2,1)
imagesc(x,z,vel);colorbar
xlabel('distance (m)');ylabel('depth (m)');
title('Exact model')
subplot(1,2,2)
imagesc(x,z,velsmo);colorbar
xlabel('distance (m)');ylabel('depth (m)');
title('Smooth mnodel')

dtstep=.0005;
itcorr=4;
tout=0:.05:tmax;
xv=x;
lap=2;
illume=1;
[shotm,shotmf,xm,zm,recsnaps,sousnaps,refsnaps,xms,zms]=rtm_shot(shot,t,xrec,zrec,xshot,zshot,w,tw,velsmo,xv,...
    dtstep,itcorr,tout,lap,illume);

r=200;
shotm=sourcemute(shotm,zm,xm,xshot,zshot,r);
shotmf=sourcemute(shotmf,zm,xm,xshot,zshot,r);

seisplottwo(shot,t,xrec,'Input shot record',shotm,zm,xm,['RTM migration (smooth model) (illume=' int2str(illume) ')']);

seisplottwo(shotm,zm,xm,'RTM migration (smooth model)',shotmf,zm,xm,'RTM migration with Laplacian filter');

plotsnaps_rtm(recsnaps,sousnaps,refsnaps,vel,xm,zm,tout,'wedgemodel','metric','jet');

%% shoot a line of 100 shots (takes a few hours) over the wedge model and then migrate the line
dx=5;
xmax=2500;zmax=1000; %maximum line length and maximum depth
vhigh=4000;vlow=2000; % high and low velocities
[vel,x,z]=wedgemodel(dx,xmax,zmax,vhigh,vlow);
fdom=35;
xshot=max(x)/3;
zshot=0.0;
dt=.004; %temporal sample rate
dtstep=.0005;%time step size
%[w,tw]=ricker(dtstep,fdom,.2);%wavelet
[w,tw]=wavemin(dt,fdom,.2);
tmax=2*zmax/vlow; %maximum time
nshots=100;
[shotsw,t,xshotsw,xrecw]=afd_shootline(dx,vel,dt,dtstep,w,tw,tmax,nshots);


%now migrate
dtstep=.0005;
itcorr=4;
tout=[];
xv=x;
lap=2;
illume=1;
shotswm=cell(size(shotsw));
shotswmf=shotswm;
[w,tw]=wavemin(dtstep,fdom,.2);
t0=clock;
for k=1:nshots
    zrec=zeros(size(xrecw{k}));
    [shotswm{k},shotswmf{k},xm,zm]=rtm_shot(shotsw{k},t,xrecw{k},zrec,xshotsw(k),0,w,tw,vel,xv,...
        dtstep,itcorr,tout,lap,illume);
    time_used=etime(clock,t0);
    time_per_shot=time_used/k;
    time_left=time_per_shot*(nshots-k);
    disp(['Finished RTM for shot ' int2str(k) ' after ' num2str(time_used/60) ' min']);
    disp(['Estimated time remaining ' num2str(time_left/60)  ' min'])
end
mute=[0,0;200,200;500,500;1000,1000];
killrad=50;
taper=10;
gainopt=0;
[stackw,gathersw,offsetsw]=migstack(shotswm,xm,zm,xshotsw,mute,killrad,gainopt);
[stackwf,gatherswf,offsetswf]=migstack(shotswmf,xm,zm,xshotsw,mute,killrad,gainopt);
seisplottwo(stackw,zm,xm,'PSRTM Wedge no filter',stackwf,zm,xm,'PSRTM Wedge Laplacian filter');

%try a butterworth low-pass to get rig of artefacts
fmin=.01;%picked by looking at shotswm{18} with seisplotfk
shotswmb=cell(size(shotswm));
for k=1:nshots
    shotswmb{k}=butterband(shotswm{k},zm,fmin,0,4,0);
end
[stackwb,gatherswb,offsetswb]=migstack(shotswmb,xm,zm,xshotsw,mute,killrad,gainopt);
seisplottwo(stackw,zm,xm,'PSRTM Wedge no filter',stackwb,zm,xm,'PSRTM Wedge Butterworth high pass filter');

stackwb2=butterband(stackw,zm,fmin,0,4,0);
seisplottwo(stackw,zm,xm,'PSRTM Wedge no filter',stackwb2,zm,xm,'PSRTM Wedge Butterworth high pass filter');

%look at CIG's
ngath=length(gathersw);
igath=round(linspace(1,ngath,25));
titles=cell(size(igath));
for k=1:length(titles)
   titles{k}=['CIG for x=' num2str(xm(igath(k)))];
end
plotgathers(gatherswf(igath),offsetswf(igath),zm,'Offset (m)','Depth (m)',titles,'PSRTM Wedge with Laplacian');

%look at some of the migrated shots
ishots=near(xshotsw,100,1500);
titles=cell(size(ishots));
offsets=titles;
for k=1:length(ishots)
   titles{k}=['PSRTM Laplacian xshot=' num2str(xshotsw(ishots(k)))];
   offsets{k}=xrecw{ishots(k)}-xshotsw(k);
end
plotgathers(shotswmf(ishots),offsets,zm,'Offset (m)','Depth (m)',titles,'Wedge migrate shots');

fnameout='wedge_rtm_shots';
save(fnameout,'shotsw','t','xshotsw','xrecw','shotswm','shotswmf','xm','zm','illume',...
    'vel','x','z','w','tw','fdom');

fnameout='wedge_rtm_stack';
save(fnameout,'stackw','stackwf','gathersw','gatherswf','offsetsw','offsetswf','t','xm','zm','illume',...
    'vel','x','z','w','tw','fdom');

fnameout='wedge_100_shots';
save(fnameout,'shotsw','t','xshotsw','xrecw','vel','x','z','w','tw','fdom');
%% shoot a line over the thrust model (50 shots) and then migrate the line. Warning, this takes a day!!
dx=5;
[velt,xt,zt]=thrustmodel(dx);
fdom=35;
dt=.004; %temporal sample rate
%[w,tw]=ricker(dtstep,fdom,.2);%wavelet
[w,tw]=wavemin(dt,fdom,.2);
tmax=2.5; %maximum time
nshots=50;
if(~exist('thrust_50_shots.mat','file'))
    dtstep=.0005;%time step size
    [shotst,tt,xshotst,xrect]=afd_shootline(dx,velt,dt,dtstep,w,tw,tmax,nshots);
    
    fnameout='Thrust_50_shots';
    save(fnameout,'shotst','tt','xshotst','xrect','velt','xt','zt','w','tw','fdom');
else
    load('Thrust_50_shots');
end

%now migrate
dtstep=.0005;
itcorr=4;
tout=[];
xv=xt;
lap=2;
illume=0;
if(illume==1)
    IClabel='ICIC';
else
    IClabel='CCIC';
end
shotstm=cell(size(shotst));
shotstmf=shotstm;
[w,tw]=wavemin(dtstep,fdom,.2);
t0=clock;
for k=1:nshots
    zrec=zeros(size(xrect{k}));
    [shotstm{k},shotstmf{k},xmt,zmt]=rtm_shot(shotst{k},tt,xrect{k},zrec,xshotst(k),0,w,tw,velt,xv,...
        dtstep,itcorr,tout,lap,illume);
    time_used=etime(clock,t0);
    time_per_shot=time_used/k;
    time_left=time_per_shot*(nshots-k);
    disp(['Finished RTM for shot ' int2str(k) ' after ' num2str(time_used/60) ' min']);
    disp(['Estimated time remaining ' num2str(time_left/60)  ' min'])
end

fnameout=['Thrust_rtm_shots_' IClabel];
save(fnameout,'shotst','tt','xshotst','xrect','shotstm','shotstmf','xmt','zmt','illume',...
    'velt','xt','zt','w','tw','fdom');

mute=[0,0;200,200;500,500;1000,1000];
killrad=50;
taper=10;
gainopt=0;
[stackt,gatherst,offsetst]=migstack(shotstm,xmt,zmt,xshotst,mute,killrad,gainopt);
[stacktf,gatherstf,offsetstf]=migstack(shotstmf,xmt,zmt,xshotst,mute,killrad,gainopt);
seisplottwo(stackt,zmt,xmt,['PSRTM Thrust no filter ' IClabel],stacktf,zmt,xmt,['PSRTM Thrust Laplacian filter ' IClabel]);

% %try a butterworth low-pass to get rig of artefacts
% fmin=.01;%picked by looking at shotswm{18} with seisplotfk
% shotstmb=cell(size(shotstm));
% for k=1:nshots
%     shotstmb{k}=butterband(shotstm{k},zmt,fmin,0,4,0);
% end
% [stacktb,gatherstb,offsetstb]=migstack(shotstmb,xmt,zmt,xshotst,mute,killrad,gainopt);
% seisplottwo(stackt,zmt,xmt,'PSRTM Thrust no filter',stacktb,zmt,xmt,'PSRTM Thrust Butterworth high pass filter');

% stacktb2=butterband(stackt,zmt,fmin,0,4,0);
% seisplottwo(stackt,zmt,xmt,'PSRTM Thrust no filter',stacktb2,zmt,xmt,'PSRTM Thrust Butterworth high pass filter');

%look at CIG's
ngath=length(gatherst);
igath=round(linspace(1,ngath,25));
titles=cell(size(igath));
for k=1:length(titles)
   titles{k}=['CIG for x=' num2str(xmt(igath(k)))];
end
plotgathers(gatherstf(igath),offsetstf(igath),zmt,'Offset (m)','Depth (m)',titles,'PSRTM Thrust with Laplacian');

%look at some of the migrated shots
ishots=near(xshotst,100,3000);
titles=cell(size(ishots));
offsets=titles;
for k=1:length(ishots)
   titles{k}=['PSRTM Laplacian xshot=' num2str(xshotst(ishots(k)))];
   offsets{k}=xrect{ishots(k)}-xshotst(k);
end
plotgathers(shotstmf(ishots),offsets,zmt,'Offset (m)','Depth (m)',titles,'Thrust migrated shots');

fnameout=['Thrust_rtm_stack_' IClabel];
save(fnameout,'stackt','stacktf','gatherst','gatherstf','offsetst','offsetstf','tt','xmt','zmt','illume',...
    'velt','xt','zt','w','tw','fdom');
%% do a cmp stack and a post stack migration of the 50 shots
fname='thrust_50_shots';
load(fname);
%make rms velocity model
dt=tt(2)-tt(1);
tmax=tt(end);
vmax=max(velt(:));
vmin=min(velt(:));
velrms=vzmod2vrmsmod(velt,zt,dt,tmax);
figure
subplot(1,2,1)
imagesc(xt,zt,velt,[vmin vmax]);colorbar
colormap(copperud)
xlabel('distance (m)');
ylabel('depth (m)')
title('Interval velocity model');
subplot(1,2,2)
imagesc(xt,tt,velrms,[vmin vmax]);colorbar
%imagesc(xt,tt,velrms);colorbar
xlabel('distance (m)');
ylabel('time (s)')
title('RMS velocity model');
prepfig

%cmp stack
dcmp=dx;
[stack,xcmp]=cmpstack(shotst,tt,xrect,xshotst,[dcmp min(xt) max(xt)],velrms,tt,xt);
load('thrust_explode_rtm');
seisplottwo(stack,tt,xcmp,'CMP stack',seisex,tt,xt,'exploding reflector model')

dtstep=.0005;
[stackm,xm2,zm2]=rtm_stack(stack,tt,xt,dtstep,velt,tout);

seisplottwo(stacktf,zm,xm,'Prestack RTM with Laplaican filter',stackm,xm2,zm2,'Post-stack RTM no filter');

seisplottwo(stacktf,zm,xm,'Prestack RTM with Laplaican filter',seism,xm2,zm2,'Exploding reflector RTM no filter');

fnameout='thrust_cmpstack_rtm';
save(fnameout,'stack','tt','xcmp','seisex','xt','zt','velt','velrms','stacktf','stackm','seism','xm','zm','xm2','zm2')
%% model and RTM migrate a single shot over the thrust model This takes 20 minutes or so
dx=5;
[velt,xt,zt]=thrustmodel(dx);
% velt=gaussian_smoother(velt,xt,zt,30);
zmax=max(zt);
fdom=35;
zshot=0.0;
%xshot=1940;
xshot=1000;
dt=.004; %temporal sample rate
dtstep=.0005;%time step size
%[w,tw]=ricker(dtstep,fdom,.2);%wavelet
[w,tw]=wavemin(dtstep,fdom,.2);
tmax=2.5; %maximum time
xrect=xt;
zrect=zeros(size(xrect));
[shott,tt]=afd_shotrec_alt(dx,dtstep,dt,tmax,velt,xshot,zshot,xrect,zrect,w,tw,2);

seisplotfk(shott,tt,xt,['Thrust model shot at xshot=' num2str(xshot)]);

% rtm_shot with exact velocity model
dtstep=.0005;
itcorr=4;
tout=0:.05:tmax;
xvt=xt;
lap=2;
illume=1;
[shottm,shottmf,xmt,zmt,recsnaps,sousnaps,refsnaps,xms,zms]=rtm_shot(shott,tt,xrect,zrect,xshot,zshot,w,tw,velt,xvt,...
    dtstep,itcorr,tout,lap,illume);

illume=0;
[shottmcc,shottmfcc,xmt,zmt,recsnapscc,sousnapscc,refsnapscc,xms,zms]=rtm_shot(shott,tt,xrect,zrect,xshot,zshot,w,tw,velt,xvt,...
    dtstep,itcorr,tout,lap,illume);

seisplottwo(shott,tt,xrect,['Thrust model shot at xshot=' num2str(xshot)],shottm,zmt,xmt,['RTM migrated shot record (illume=' int2str(illume) ')']);

seisplottwo(shottmf,zmt,xmt,'RTM migration with Laplacian filter',shottm,zmt,xmt,'RTM migrated shot record');

seisplottwo(shottm,zmt,xmt,'RTM migrated shot, DECIC',shottmcc,zmt,xmt,'RTM migrated shot, CCIC');
seisplottwo(shottmf,zmt,xmt,'RTM migrated shot, Laplacian filter, DECIC',shottmfcc,zmt,xmt,'RTM migrated shot, Laplacian filter, CCIC');

% fmin=.01;%picked by looking at shot with seisplotfk
% shottmb=butterband(shottm,zmt,fmin,0,4,0);
% seisplottwo(shottmb,zmt,xmt,'RTM migration with butterband highpass',shottmf,zmt,xmt,'RTM migration with Laplacian filter');


plotsnaps_rtm(recsnaps,sousnaps,refsnaps,velt,xmt,zmt,tout,'Thrust model','metric','copperud')

plotsnaps_rtm_big(recsnaps,sousnaps,refsnaps,velt,nan,nan,tout,'Thrust model','metric','copperud',[1 1 4])

fnameout=['RTM_thrust_snapshots_x=' int2str(xshot)];
save(fnameout,'shott','tt','xt','xshot','velt','zt','xrect','zrect','w','tw','fdom','shottm',...
    'shottmf','xmt','zmt','recsnaps','sousnaps','refsnaps','xms','zms','tout')