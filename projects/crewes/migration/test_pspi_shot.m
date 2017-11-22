%% model and RTM migrate a single shot over the thrust model This thaks 20 minutes or so
dx=5;
[velt,xt,zt]=thrustmodel(dx);
% velt=gaussian_smoother(velt,xt,zt,30);
zmax=max(zt);
fdom=35;
zshot=0.0;
xshot=4000;%change this in both cells 1&2 to move the shot
dt=.004; %temporal sample rate
dtstep=.0005;%time step size
%[w,tw]=ricker(dtstep,fdom,.2);%wavelet
[w,tw]=wavemin(dtstep,fdom,.2);
tmax=2.5; %maximum time
xrect=xt;
zrect=zeros(size(xrect));
[shott,tt]=afd_shotrec_alt(dx,dtstep,dt,tmax,velt,xshot,zshot,xrect,zrect,w,tw,2);

figure
imagesc(xt,zt,velt);colorbar,colormap(copperud)
xlabel('distance (m)')
ylabel('depth (m)')
title('Thrust model, * marks the shot')
line(xshot,0,'color','r','linestyle','none','marker','*','markersize',12)
prepfig

seisplotfk(shott,tt,xrect,['Thrust model shot at xshot=' num2str(xshot)]);

save(['thrust_shot_x=' int2str(xshot)],'shott','tt','xrect','zrect','velt','xt','zt','fdom');

%% migrate the shot record 
%this if block looks for the saved dataset from the previous cell
xshot=4000;
if(~exist('seistt','var'))
    if(exist(['thrust_shot_x=' int2str(xshot) '.mat'],'file'))
        load(['thrust_shot_x=' int2str(xshot)])
    else
        disp('Please run the first cell of this script')
    end
end


frange=[0 100];
stab=.001;
zout=0:50:2100;
[w,tw]=wavemin(dt,fdom,.2);
[shotmigdec,shotmigcc,illumination,chkpnts]=pspi_shot(shott,tt,xrect,velt,xt,zt,xshot,frange,stab,zout,w,tw);

seisplottwo(shott,tt,xrect,'Thrust model unmigrated shot',shotmigcc,zt,xt,'Migrated shot CCIC');

seisplottwo(shotmigdec,zt,xrect,'Migrated shot DECIC',shotmigcc,zt,xt,'Migrated shot CCIC');

%plotchkpnts_pspi(chkpnts,tt,velt,xt,zt,'Thrust model','metric','jet');

plotchkpnts_pspi_big(chkpnts,velt,'Thrust model','copperud',[1 1 8 1],4);

%% model and migrate a 50 shot test line
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
frange=[0 100];
stab=.001;
shotsdec=cell(size(shotst));
shotscc=shotsdec;
t0=clock;
for k=1:nshots
    [shotsdec{k},shotscc{k}]=pspi_shot(shotst{k},tt,xrect{k},velt,xt,zt,xshotst(k),frange,stab);
    time_used=etime(clock,t0);
    time_per_shot=time_used/k;
    time_left=time_per_shot*(nshots-k);
    disp(['Finished RTM for shot ' int2str(k) ' after ' num2str(time_used/60) ' min']);
    disp(['Estimated time remaining ' num2str(time_left/60)  ' min'])
end

fnameout='Thrust_pspi_shots';
save(fnameout,'shotst','tt','xshotst','xrect','shotsdec','shotscc',...
    'velt','xt','zt','frange','stab');

mute=[0,0;200,200;500,500;1000,1000];
killrad=50;
taper=10;
gainopt=0;
[stackdec,gathersdec,offsetsdec]=migstack(shotsdec,xt,zt,xshotst,mute,killrad,gainopt);
[stackcc,gatherscc,offsetscc]=migstack(shotscc,xmt,zmt,xshotst,mute,killrad,gainopt);
seisplottwo(stackdec,zmt,xmt,'PSPI Thrust DECIC',stackcc,zmt,xmt,'PSPI Thrust CCIC');


% %look at CIG's
% ngath=length(gatherst);
% igath=round(linspace(1,ngath,25));
% titles=cell(size(igath));
% for k=1:length(titles)
%    titles{k}=['CIG for x=' num2str(xmt(igath(k)))];
% end
% plotgathers(gatherstf(igath),offsetstf(igath),zmt,'Offset (m)','Depth (m)',titles,'PSRTM Thrust with Laplacian');
% 
% %look at some of the migrated shots
% ishots=near(xshotst,100,3000);
% titles=cell(size(ishots));
% offsets=titles;
% for k=1:length(ishots)
%    titles{k}=['PSRTM Laplacian xshot=' num2str(xshotst(ishots(k)))];
%    offsets{k}=xrect{ishots(k)}-xshotst(k);
% end
% plotgathers(shotstmf(ishots),offsets,zmt,'Offset (m)','Depth (m)',titles,'Thrust migrated shots');

fnameout='Thrust_prestack_pspi_stack';
save(fnameout,'stackdec','stackcc','gathersdec','gatherscc','offsetsdec','offsetscc','tt',...
    'velt','xt','zt');

%% load and examine the stacks
fname='Thrust_prestack_pspi_stack.mat';
load(fname)

seisplottwo(stackdec,zt,xt,'PSPI Thrust DECIC',stackcc,zt,xt,'PSPI Thrust CCIC');