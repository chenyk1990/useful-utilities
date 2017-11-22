%demo ray migration on a series of models
%% #1 this cell creates a single diffraction
%constant velocity
modelnamed='A single diffraction';
xmax=2000;
zmax=1000;
xnot=xmax/2;%x coordinate of the diffractor
znot=500;%z coordinate of the diffractor
v=2000;
dx=5;
dt=.002;
tmax=zmax*2/v;
xd=0:dx:xmax;
td=0:dt:tmax;
zd=0:dx:zmax;
filtparms=[5,10,80,100]; %filter specification Ormsby style
%Create the diffraction
tnot=2*znot/v;
seis=zeros(length(td),length(xd));
seis=event_hyp(seis,td,xd,tnot,xnot,v,1);
seisdiff=filtf(seis,td,[filtparms(2) filtparms(2)-filtparms(1)],[filtparms(3) filtparms(4)-filtparms(3)]);
% create the depth section
depth=zeros(length(zd),length(xd));
ix=near(xd,xnot);
iz=near(zd,znot);
depth(iz,ix)=1;
g=gaus_radial(dx,5*dx);
depth=conv2(depth,g,'same');
veldiff=v*ones(size(depth));
Amax=max(depth(:));
ind=find(depth>.5*Amax);
veldiff(ind)=1.5*v;

raymig(seisdiff,veldiff,td,xd,zd,modelnamed)

%% #2 Very simple dipping reflector
%dipping reflector
%do a finite-difference model of dipping reflector
modelnamedip='dip model';
dx=5;
vlow=2000;vhigh=4000;
xmax=2500;zmax=1200;
dip=10;
[veldip,xdip,zdip]=dipmodel(dx,xmax,zmax,vhigh,vlow,dip);

%Make a ZOS
%exploding reflector
dt=.004; %temporal sample rate
dtstep=.001;
tmax=2*zmax/vlow; %maximum time
[seisfiltd,seis,tdip]=afd_explode(dx,dtstep,dt,tmax, ...
 		veldip,xdip,zeros(size(xdip)),[5 10 40 50],0,2);

raymig(seisfiltd,veldip,tdip,xdip,zdip,modelnamedip)

%% #3 Ray migration of the buried cosyncline (interactive)
%first remake the cosyncline with ivery=1;
ievery=1;%make a point diffraction every this many traces
polarity=1;%use a 1 to get a major syncline in the middle and anticlines on the flanks. A -1 gets the opposite
znot=500;%average depth of the cosine, try values like 500 and 150
amp=200;%amplitude of cosine

xperiod=1000;
xmax=2000;
tmax=1;
xnot=500;%x coordinate of zero argument


v=2000;
dx=5;
dt=.002;
x=0:dx:xmax;
t=0:dt:tmax;

filtparms=[5,10,80,100]; %filter specification Ormsby style

%diffraction coordinates
xd=x(1:ievery:end);
zd=znot-polarity*amp*cos((xd-xnot)*2*pi/xperiod);
td=2*zd/v;

%ok install diffractions
seis=zeros(length(t),length(x));
for k=1:length(xd)
    seis=event_hyp(seis,t,x,td(k),xd(k),v,1);
end

%filter
seis=filtf(seis,t,[filtparms(2) filtparms(2)-filtparms(1)],[filtparms(3) filtparms(4)-filtparms(3)]);

%do a ray migration of the previous
zmax=1000;
z=0:dx:zmax;
zd=td*v/2;
vel=v*ones(length(z),length(x));
for k=1:length(x)
    ind=find(z>zd(k));
    vel(ind,k)=3000;
end

[h1,h2,h3]=raymig(seis,vel,t,x,z,'Cosyncline');


%% #4 this cell uses the wedge model. It is useful to understand depth migration
% Try migrating pickc from the top of the anticline with both the exact and smooth models.
% Which model give better results? Why?
%

%do a finite-difference model of wedge
modelnamew='wedge model';
dx=5;
vlow=2000;vhigh=4000;
xmax=2500;zmax=1200;
[velw,xw,zw]=wedgemodel(dx,xmax,zmax,vhigh,vlow);
dt=.004; %temporal sample rate
dtstep=.001;
tmax=2*zmax/vlow; %maximum time
[seisfiltw,seis,tw]=afd_explode(dx,dtstep,dt,tmax, ...
 		velw,xw,zeros(size(xw)),[5 10 40 50],0,2);
    
raymig(seisfiltw,velw,tw,xw,zw,modelnamew)


%% #5 THis cell uses the channel model. It is boring except for the channel diffraction
%do a finite-difference model of channel
modelnamec='channel model';
dx=5;
vlow=2000;vhigh=4000;
xmax=2500;zmax=1200;
[velc,xc,zc]=channelmodel(dx,xmax,zmax,vhigh,vlow);
dt=.004; %temporal sample rate
dtstep=.001;
tmax=2*zmax/vlow; %maximum time
[seisfiltc,seis,tc]=afd_explode(dx,dtstep,dt,tmax, ...
 		velc,xc,zeros(size(xc)),[5 10 40 50],0,2);
    
raymig(seisfiltc,velc,tc,xc,zc,modelnamec)

%% #6 The thrust model has a shadow zone beneath the thrust
%do a finite-difference model of thrust
modelnamet='thrust model';
dx=5;
vlow=2000;vhigh=3500;
xmax=5100;zmax=2500;
[velt,xt,zt]=thrustmodel(dx,xmax,zmax,vhigh,vlow);
dt=.004; %temporal sample rate
dtstep=.001;
tmax=2*zmax/vlow; %maximum time
[seisfiltt,seis,tt]=afd_explode(dx,dtstep,dt,tmax, ...
 		velt,xt,zeros(size(xt)),[5 10 40 50],0,2);
    
raymig(seisfiltt,velt,tt,xt,zt,modelnamet)


%% #7 depth migration of thrust using PSPI
% this is nice to compare to the raytracing in the previous cell

usesmooth=0;%0 to migrate with the exact model, 1 to use the smoothed model

if(usesmooth==1)
    zcheck=0:100:2000;
    % [zosmig,exzos]=pspi_stack(seisfiltt,t,x,velt,x,z,[5 50],zcheck);
    [zosmigsmo,exzossmo]=pspi_stack(seisfiltt,t,x,veltsmo,x,z,[5 50],zcheck);

    plotimage(zosmigsmo,x,z);
    title('Migrated with smoothed model')
    
    xs=cell(size(exzossmo));
    ts=cell(size(exzossmo));
    titles=cell(size(exzossmo));
    for k=1:length(exzossmo)
        xs{k}=(x(2)-x(1))*(0:size(exzossmo{k},2)-1);
        ts{k}=(t(2)-t(1))*(0:size(exzossmo{k},1)-1);
        titles{k}=['Extrapolated (smooth model) to depth ' int2str(zcheck(k)) 'm'];
    end
    
    %load the extrapolations into plotgathers
    plotgathers(exzossmo,xs,ts,'distance (m)','time (s)',titles);
else
    zcheck=0:100:2000;
    [zosmig,exzos]=pspi_stack(seisfiltt,t,x,velt,x,z,[5 50],zcheck);
    
    plotimage(zosmig,x,z);
    title('Migrated with exact model')
    
    xs=cell(size(exzos));
    ts=cell(size(exzos));
    titles=cell(size(exzos));
    for k=1:length(exzos)
        xs{k}=(x(2)-x(1))*(0:size(exzos{k},2)-1);
        ts{k}=(t(2)-t(1))*(0:size(exzos{k},1)-1);
        titles{k}=['Extrapolated to depth ' int2str(zcheck(k)) 'm'];
    end

    %load the extrapolations into plotgathers
    plotgathers(exzos,xs,ts,'distance (m)','time (s)',titles);
end

