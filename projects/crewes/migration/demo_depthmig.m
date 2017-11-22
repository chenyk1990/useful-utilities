%A series of depth versus time migration examples
% Sections 1-5 use the Thrust model to compare time and depth migration. This is quite a complex model
% Sections 6-9 use the much simpler wedge model in a similar comparison.

%% #1 do a finite-difference model of thrust

modelname='thrust model';
dx=5;
vlow=2000;vhigh=3500;
xmax=5100;zmax=2500;
[velt,x,zv]=thrustmodel(dx,xmax,zmax,vhigh,vlow);
figure;imagesc(x,zv,velt);colorbar
pititle('The thrust model')
dt=.004; %temporal sample rate
dtstep=.001;
tmax=2*zmax/vlow; %maximum time
[seisfiltt,seis,t]=afd_explode(dx,dtstep,dt,tmax, ...
 		velt,x,zeros(size(x)),[5 10 40 50],0,2);
    
seisplot(seisfiltt,t,x,'Thrust model unmigrated');

%% #2 depth migration of thrust

zcheck=0:100:2000;
[zosmig_depth,exzos]=pspi_stack(seisfiltt,t,x,velt,x,zv,[5 50],zcheck);
seisplot(zosmig_depth,zv,x,'Depth migration of Thrust');
xs=cell(size(exzos));
ts=cell(size(exzos));
titles=cell(size(exzos));
for k=1:length(exzos)
    xs{k}=(x(2)-x(1))*(0:size(exzos{k},2)-1);
    ts{k}=(t(2)-t(1))*(0:size(exzos{k},1)-1);
    titles{k}=['Extrapolated to ' int2str(zcheck(k))];
end
save thrustdata

%load the extrapolations into plotgathers
plotgathers(exzos,xs,ts,'distance (m)','time (s)',titles);

%% #3 Kirchoff time migration of thrust
%build an RMS velocity model
[vrms,tv]=vzmod2vrmsmod(velt,zv,dt,tmax);

%use kirk because its faster than kirk_mig
params=nan*(ones(1,5));
[zosmig_timek,tmig,xmig]=kirk_mig(seisfiltt,vrms,t,x,params);

seisplot(zosmig_timek,tmig,xmig,'Kirchhoff time migration of thrust')

%% #4 PSPI time migration of thrust
%build an RMS velocity model
[vrms,tv]=vzmod2vrmsmod(velt,zv,dt,tmax);

%use kirk because its faster than kirk_mig
params=nan*(ones(1,5));
taucheck=0:.1:2.;
[zosmig_timep,ex_zos]=pspi_stack_tmig_rms(seisfiltt,t,x,vrms,x,tv,[5 50],taucheck);

seisplot(zosmig_timep,tv,x,'PSPI time migration of thrust')

xs=cell(size(exzos));
ts=cell(size(exzos));
titles=cell(size(exzos));
for k=1:length(exzos)
    xs{k}=(x(2)-x(1))*(0:size(exzos{k},2)-1);
    ts{k}=(t(2)-t(1))*(0:size(exzos{k},1)-1);
    titles{k}=['Extrapolated to ' int2str(taucheck(k))];
end

%load the extrapolations into plotgathers
plotgathers(exzos,xs,ts,'distance (m)','time (s)',titles);

%% #5 raytrace migration of thrust

raymig(seisfiltt,velt,t,x,zv,'Thrust model')

%% #6 do a finite-difference model of wedge
modelnamew='Wedge model';
dx=2;
vlow=2000;vhigh=4000;
xmax=2500;zmax=1200;
[velw,xw,zw]=wedgemodel(dx,xmax,zmax,vhigh,vlow);
dt=.004; %temporal sample rate
dtstep=.0005;
tmax=2*zmax/vlow; %maximum time
[seisfiltw,seis,tw]=afd_explode(dx,dtstep,dt,tmax, ...
    velw,xw,zeros(size(xw)),[5 10 40 50],0,2);

[w,tww]=ricker(dtstep,30,.2);
%[w,tww]=ormsby(5,10,40,50,.2,dtstep);
[seisfiltwa,t,r]=afd_explode_alt(dx,dtstep,dt,tmax, ...
    velw,xw,zeros(size(xw)),w,tww,2);

seisplot(seisfiltwa,t,xw,modelnamew);
figure;imagesc(xw,zw,velw);colorbar
pititle(modelnamew);

%% #7 depth migration 
zosmig_depth=pspi_stack(seisfiltw,tw,xw,velw,xw,zw,[5 50]);
seisplot(zosmig_depth,zw,xw,['Depth migration of ' modelnamew]);



%% #8 Kirchhoff time migration
%build an RMS velocity model
[vrmsw,tvw]=vzmod2vrmsmod(velw,zw,dt,tmax);
%use kirk because its faster than kirk_mig
params=nan*(ones(1,5));
[zosmig_timekw,tmigw,xmigw]=kirk_mig(seisfiltw,vrmsw,tw,xw,params);

seisplot(zosmig_timekw,tmigw,xmigw,['Kirchhoff time migration of ' modelnamew])

%% #9 raytrace migration of wedge
raymig(seisfiltw,velw,tw,xw,zw,modelnamew)