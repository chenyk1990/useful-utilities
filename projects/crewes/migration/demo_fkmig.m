% a series of demonstrations of f-k migration on various models
%%
% #1 migrate a single diffraction
%
% try changing fmax and dip max and observe the consequences
%
% Since migration is a linear process, if we consider a seismic section to be the sum of many
% many diffractions, then this is really the whole story. 
% How does dipmax affect the size of the focal point? How about fmax? Does changing fmax only
% change the vertical size or does the lateral size change also? Same question for dipmax? If
% the goal is optimal resolution, then what values do we want for fmax and dipmax?
% Try adding random noise to the data (set s2n to a number like 2) and see how this affects
% things.
%



fmax=50; %maximum sigma frequency
dipmaxs=[80 30];%maximum dip to migrate
%s2n=.5;%signal to noise ratio. 
s2n=inf;%this means infinite signal to noise and is how you turn the noise to zero
dt=.004;
tmax=2.0;
dx=5;
nx=1000;
v=2000;
x=(0:nx-1)*dx;
nt=round(tmax/dt)+1;
t=(0:nt-1)*dt;
seis=zeros(nt,nx);
tnot=.5;
xnot=max(x)/2;
seis=event_hyp(seis,t,x,tnot,xnot,v,1);
%filter
seis=filtf(seis,t,[10 5],[fmax 10],0);
if(~isinf(s2n))
    noise=rnoise(seis(:),s2n);
    seis=seis+reshape(noise,nt,nx);
end
params=nan*ones(1,13);
params(1)=fmax;%frequency limit
params(3)=dipmaxs(1);%dip limit
[seismig1,tmig,xmig]=fkmig(seis,t,x,v,params);
params(3)=dipmaxs(2);%dip limit
[seismig2,tmig,xmig]=fkmig(seis,t,x,v,params);
zz=tmig*v/2;

data=seisplotfk(seis,t,x,'Diffraction hyperbola');
axes(data{2});
f=data{3};
k=data{4};
fbdy=v*abs(k)/2;
ind=find(fbdy<=max(f));
h=line(k(ind),fbdy(ind),'color','r');
legend(h,'evanescent boundary')


data2=seisplotfk(seismig1,zz,x,['f-k migration of diffraction, dipmax= ' int2str(dipmaxs(1)) 'deg']);
axes(data2{2});
yl=get(gca,'ylim');
axis equal
ylim(yl)

data2=seisplotfk(seismig2,zz,x,['f-k migration of diffraction, dipmax= ' int2str(dipmaxs(2)) 'deg']);
axes(data2{2});
yl=get(gca,'ylim');
axis equal
ylim(yl)


%expanded views of focal points
delt=.1;
delx=100;
indt=near(t,tnot-delt,tnot+delt);
indx=near(x,xnot-delx,xnot+delx);
seisplottwo(seismig1(indt,indx),t(indt),x(indx),['Migrated diffraction, dipmax= ' int2str(dipmaxs(1)) 'deg'],...
    seismig2(indt,indx),t(indt),x(indx),['Migrated diffraction, dipmax= ' int2str(dipmaxs(2)) 'deg']);
%% #2
% use fkmig to show the impulse response of migration
% try changing fmax and dip max and observe the consequences
%
% Why is this called the impulse response of migration? What shape is the impulse response?
% does the shape of the response change with the position of the impulse?
%
% Why does the unmigrated spectrum extend outside the evanescent boundaries? Is this
% meaningful?

fmax=75;%maximum frequency to me migrated
dipmax=90;%maximum dip to be migrated
tmax=2.0;%maximum time
xmax=2500;%line length
tnot=.75*tmax;%time of the impulse
xnot=.5*xmax;%x location of the impulse
dt=.004;

dx=10;
nx=round(xmax/dx)+1;
v=2000;
x=(0:nx-1)*dx;
nt=round(tmax/dt)+1;
t=dt*(0:nt-1)';
seis=zeros(nt,nx);
it=near(t,tnot);
ix=near(x,xnot);
tmp=impulse(t,it(1));
tmp=filtf(tmp,t,[10 5],[50 10]);
seis(:,ix(1))=tmp;
params=nan*ones(1,13);
params(1)=fmax;%maximum frequency to migrate
params(3)=dipmax;%maximum dip to migrate
[seismig,tmig,xmig]=fkmig(seis,t,x,v,params);
zmig=v*tmig/2;

data1=seisplotfk(seis,t,x,'Impulse to be migrated');
axes(data1{1})
h=line(x(ix),t(it),'linestyle','none','marker','o','color','k','markersize',12);
legend(h,'Impulse location')
axes(data1{2})
f=data1{3};
kx=data1{4};
%draw evanescent boundaries
fbdy=v*abs(kx)/2;
ind=find(fbdy<=max(f));
h=line(kx(ind),fbdy(ind),'color','r');
legend(h,'evanescent boundary')

data2=seisplotfk(seismig,zmig,xmig,['Migrated impulse, Dipmax=' num2str(dipmax)]);

axes(data2{2});
yl=get(gca,'ylim');
axis equal
%set(gca,'ylim',yl,'xlim',[-.05 .05]);
set(gca,'ylim',yl);

%% #3
% demo migration of a single live trace
% try changing fmax and dip max and observe the consequences
%
% The migration of a section of traces is just the superposition of many images like this, one
% for each trace (because migration is a linear process). The lateral extent of this respnse
% gives an idea of how far data can move during migration.
%
fmax=75;
dipmax=80;
dt=.004;
tmax=2.0;
xmax=5000;
nt=round(tmax/dt)+1;
dx=10;
nx=round(xmax/dx)+1;
nz=500;
v=2000;
x=(0:nx-1)*dx;
seis=zeros(nt,nx);
[r,t]=reflec(tmax,dt);

s=filtf(r,t,[10 5],[50 10]);
% w=ricker(dt,40);
% s=convz(r,w);
seis(:,round(nx/2))=s(1:nt);
params=nan*ones(1,13);
params(1)=fmax;
params(3)=dipmax;%specify dip limit
[seismig,tmig,xmig]=fkmig(seis,t,x,v,params);
zmig=v*tmig/2;

data1=seisplotfk(seis,t,x,'Trace to be migrated');

axes(data1{2})
f=data1{3};
kx=data1{4};
%draw evanescent boundaries
fbdy=v*abs(kx)/2;
ind=find(fbdy<=max(f));
h=line(kx(ind),fbdy(ind),'color','r');
legend(h,'evanescent boundary')

data2=seisplotfk(seismig,zmig,xmig,['Migrated trace, dipmax=' num2str(dipmax)]);
axes(data2{1});
yl=get(gca,'ylim');
axis equal
%set(gca,'ylim',yl,'xlim',[-.05 .05]);
set(gca,'ylim',yl);
axes(data2{2});
yl=get(gca,'ylim');
axis equal
%set(gca,'ylim',yl,'xlim',[-.05 .05]);
set(gca,'ylim',yl);

%% #4
% make a standard synthetic to migrate
% try changing fmax and dip max and observe the consequences
%
% How does cahnging dipmax from 80 to 30 change the result? Is the only effect the loss of
% steep disp or is there a loss of resolution as well? Whay?
%

fmax=75;
dipmax=80;
dt=.002;
dx=5;
tmax=2;
xmax=4000;
v=2500;
[w,tw]=ricker(dt,40,.2);
[seis,t,x]=makestdsynh(dt,dx,tmax,xmax,v,w,tw);
%fkmig
params=nan*ones(1,13);
params(1)=fmax;
params(3)=dipmax;
[seismig,tmig,xmig]=fkmig(seis,t,x,v,params);
zmig=v*tmig/2;

%plot input
data1=seisplotfk(seis,t,x,'Standard synthetic displayed in time');
axes(data1{1})
h1=line(max(x)/3,max(t)*.167,'color','r','linestyle','none','marker','o');
line(2*max(x)/3,max(t)*.75,'color','r','linestyle','none','marker','o');
legend(h1,'noise spikes')
axes(data1{2})
f=data1{3};
k=data1{4};
%draw evanescent boundaries
fbdy=v*abs(k)/2;
ind=find(fbdy<=max(f));
h=line(k(ind),fbdy(ind),'color','r');
legend(h,'evanescent boundary')

%Plot output
data2=seisplotfk(seismig,zmig,x,['Migrated standard synthetic in depth, dipmax=' num2str(dipmax)]);

%seisplot(seis,t,x,{'Standard section with default parameters','

%% #5 Illustrate the importance of modelling with diffractions
dt=.002;%time sample size
dx=5;%space sample size
xmax=2000;%line length
tmax=1;%record length
v=2000;%velocity
znot=300;%reflector depth
tnot=2*znot/v;%reflector time
t=(0:dt:tmax)';
x=0:dx:xmax;
%make a reflector with holes using diffractions
seish=zeros(length(t),length(x));seis=seish;
seish=event_diph(seish,t,x,v,0,500,znot,0,1);
seish=event_diph(seish,t,x,v,600,800,znot,0,1);
seish=event_diph(seish,t,x,v,850,1000,znot,0,1);
seish=event_diph(seish,t,x,v,1030,1150,znot,0,1);
seish=event_diph(seish,t,x,v,1160,1260,znot,0,1);
seish=event_diph(seish,t,x,v,1265,2000,znot,0,1);
seish=filtf(seish,t,[10 5],[100 20]);
%Make the same model but without diffractions
seis=event_dip(seis,t,x,[tnot tnot],[0 500]);
seis=event_dip(seis,t,x,[tnot tnot],[600 800]);
seis=event_dip(seis,t,x,[tnot tnot],[850 1000]);
seis=event_dip(seis,t,x,[tnot tnot],[1030 1150]);
seis=event_dip(seis,t,x,[tnot tnot],[1160 1260]);
seis=event_dip(seis,t,x,[tnot tnot],[1265 2000]);
seis=filtf(seis,t,[10 5],[100 20]);

seisplottwo(seish,t,x,'Reflection with holes and diffractions',seis,t,x,'Reflector with holes but no diffractions');

%fkmigrate both
seishm=fkmig(seish,t,x,v);
seism=fkmig(seis,t,x,v);

seisplottwo(seishm,t,x,'Migration with diffractions',seism,t,x,'Migration without diffractions');

%extract amplitudes

ind=near(t,tnot);
ind2=near(x,1600);
ah=seishm(ind(1),:);
a=seism(ind(1),:);
figure
plot(x,ah/ah(ind2),x,a/a(ind2));
xlabel('distance (m)')
ylabel('amplitude')
ylim([0 1.2])
legend('Migration with diffractions','Migration without diffractions')
title('Extracted amplitudes from migration with and without diffractions');
grid
prepfiga

%zooms around tiny hole
indt=near(t,tnot-delt,tnot+delt);
indx=near(x,1110,1200);

seisplottwo(seishm(indt,indx),t(indt),x(indx),'Tiny hole with diffractions',...
    seism(indt,indx),t(indt),x(indx),'Tiny hole without diffractions');

%% #6 Diffraction section  constant velocity simulation
%make a section full of diffractions, 
% try changing fmax and dip max and observe the consequences
%
% Study how resolution changes as you move around the section. Where is resolution the best?
% The worst? What always happens to resolution at the bottom oaf a section?
%
dx=10;
dt=.004;
xmax=2000;
zmax=2000;%tmax is computed automatically as the 2-way vertical traveltime to zmax
fmax=75;%max frequency to care about
fmin=10;%min frequency to care about
v0=2000;
c=0;

[seis,t,x]=diffraction_section(dx,dt,xmax,zmax,v0,c,fmin,fmax);
dipmax=80;
fmin2=[10 5];fmax2=[fmax 10];
nt=length(t);nx=length(x);

%fkmig
params=nan*ones(1,13);
params(1)=fmax;
params(3)=dipmax;
[seismig,tmig,xmig]=fkmig(seis,t,x,v0,params);

data1=seisplotfk(seis,t,x,['Diffraction section, v0=' num2str(v0) ' constant']);
axes(data1{2});
f=data1{3};
k=data1{4};
%draw evanescent boundaries
fbdy=v0*abs(k)/2;
ind=find(fbdy<=max(f));
h=line(k(ind),fbdy(ind),'color','r');
legend(h,'evanescent boundary')

data2=seisplotfk(seismig,tmig,xmig,['Migrated diffraction section , v0=' num2str(v0) ' constant']);

%% #7 Diffraction section  v(z) simulation
%make a section full of diffractions, 
% try changing fmax and dip max and observe the consequences
%
% Compare this to the constant velocity case and note any differences. Which is more realistic?
%
dx=10;
dt=.004;
xmax=2000;
zmax=2000;%tmax is computed automatically as the 2-way vertical traveltime to zmax
fmax=75;%max frequency to care about
fmin=10;%min frequency to care about
v0=2000;
c=.6;

[seis,t,x,v]=diffraction_section(dx,dt,xmax,zmax,v0,c,fmin,fmax);
dipmax=80;
fmin2=[10 5];fmax2=[fmax 10];
nt=length(t);nx=length(x);

%fkmig
params=nan*ones(1,13);
params(1)=fmax;
params(3)=dipmax;
[seismig,tmig,xmig]=vz_fkmig(seis,t,x,v,params);

data1=seisplotfk(seis,t,x,['Diffraction section, v0=' num2str(v0) ', c=' num2str(c)]);
axes(data1{2});
f=data1{3};
k=data1{4};
%draw evanescent boundaries
fbdy=v0*abs(k)/2;
ind=find(fbdy<=max(f));
h=line(k(ind),fbdy(ind),'color','r');
legend(h,'evanescent boundary')

data2=seisplotfk(seismig,tmig,xmig,'Migrated diffraction section v(z)');
%% #8 buried syncline. This will be a circular depression imposed on a horizontal reflector
%
%                                   0
%
% ---------------------------               -----------------------------
%                             *           *
%                               **     **
%                                 *****
%
% Notice how ievery controls the width of the unmigrated spectrum. With
% ievery=1, the unmigrated spectrum does not seem to extend to the
% evanescent boundaryies. But with ievery=5, it does. Can you explain this?
%
ievery=1;%make a point diffraction every this many traces
fmax=75;%maximum frequency to migrate
dipmax=80;%maximum fip to migrate
v=2000; %velocity of simulation
dt=.004;%time sample rate
dx=5;%x (cmp) sample size
tmax=1;%maximum time
xmax=2000;%maximum x (linelength)
zmax=1000;
x=0:dx:xmax;
t=0:dt:tmax;
z=t*v/2;
dz=z(2);

filtparms=[5,10,fmax-10,fmax]; %filter specification Ormsby style


zr=zmax/3;%depth to reflector

%x0,z0 are coordinates of the focus (the 0 in the diagram)
x0=xmax/2;
%try zr/2 and -zr/2 for the following.
z0 = zr/2;%Should be no larger than zr. A positive number for z0 means a 
%buried syncline and the seismic image will have a reverse-time branch. A
%negative number means the syncline is not buried and there is no reverse
%time branch

r=1.5*(zr-z0);%radius of the syncline. This should be greater than zr-z0 or there will be no syncline at all
x1=x0-sqrt(r^2-(zr-z0)^2);%start of syncline
x2=x0+sqrt(r^2-(zr-z0)^2);%end of syncline

%diffraction coordinates
xd=x(1:ievery:end);
zd=zeros(size(xd));
ind=xd<=x1;
zd(ind)=zr;
ind=xd>=x2;
zd(ind)=zr;
ind=find(zd==0);
for k=1:length(ind)
    xdtmp=xd(ind(k));
    zd(ind(k))=z0+sqrt(r^2-(x0-xdtmp)^2);
end
td=2*zd/v;

%ok install diffractions
seis=zeros(length(t),length(x));
for k=1:length(xd)
    seis=event_hyp(seis,t,x,td(k),xd(k),v,1);
end

%filter
seis=filtf(seis,t,[filtparms(2) filtparms(2)-filtparms(1)],[filtparms(3) filtparms(4)-filtparms(3)]);



%fkmig
params=nan*ones(1,13);
params(1)=fmax;
params(3)=dipmax;
[seismig,tmig,xmig]=fkmig(seis,t,x,v,params);
zmig=v*t/2;

%plot input
data1=seisplotfk(seis,t,x,'Buried syncline',fmax);
axes(data1{2});
f=data1{3};k=data1{4};
%draw evanescent boundaries
fbdy=v*abs(k)/2;
ind=find(fbdy<=max(f));
h=line(k(ind),fbdy(ind),'color','r');
legend(h,'evanescent boundary')

%plot output
data2=seisplotfk(seismig,zmig,x,'Migrated buried syncline',2*fmax/v);
axes(data2{1})
h=line(xd,.5*v*td,'color','r');
set(h,'linestyle',':','linewidth',2)
h2=line(x0,z0,'color','r');
set(h2,'linestyle','none','marker','o')
legend([h h2],'Actual syncline','focal point')
axes(data2{2})
yl=get(gca,'ylim');


%% #9 model a buried cosine or "cosyncline"
%
% Notice the the unmigrated data only extends about halfway to the evanescent boundary. Why do
% you think this is so"
%
% Also notice that parts of the flanks of the cosyncline are not imaged after migration. Can
% you explain this?
%
ievery=1;%put a diffraction every this many traces
polarity=1;%use a 1 to get a major syncline in the middle and anticlines on the flanks. A -1 gets the opposite
znot=500;%mean depth of the cosine, try values like 500 and 150
fmax=75;%maximum signal frequency
dipmax=80;%maximum dip

xperiod=1000;
xmax=2000;
tmax=1;
xnot=500;%x coordinate of zero argument
amp=100;%amplitude of cosine
% zmax=1000;
v=2000;
dx=5;
dt=.004;
x=0:dx:xmax;
t=0:dt:tmax;

filtparms=[5,10,fmax-10,fmax]; %filter specification Ormsby style

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

%fkmig
params=nan*ones(1,13);
params(1)=fmax;
params(3)=dipmax;
[seismig,tmig,xmig]=fkmig(seis,t,x,v,params);
zmig=v*tmig/2;

%plot the input
data1=seisplotfk(seis,t,x,'Cosyncline unmigrated');
axes(data1{2});
f=data1{3};
k=data1{4};
%draw evanescent boundaries
fbdy=v*abs(k)/2;
ind=find(fbdy<=max(f));
h=line(k(ind),fbdy(ind),'color','r');
legend(h,'evanescent boundary')

%plot the output
data2=seisplotfk(seismig,zmig,xmig,'Migrated cosyncline');
axes(data2{1})
h=line(xd,zd,'color','r');
set(h,'linestyle',':','linewidth',2)
legend(h,'Actual cosyncline')


%% #10 dipping reflector
%
% Notice how ievery controls the width of the unmigrated spectrum. With
% ievery=1, the unmigrated spectrum does not seem to extend to the
% evanescent boundaryies. But with ievery=5, it does. Can you explain this?
%
clear
fmax=75;
dipmax=80;
ievery=1;
xmax=2000;
zmax=1000;
xnot=1500;%x coordinate where the reflector outcrops
z1=100;%first depth with nonzero amplitude
z2=500;%last depth with nonzero amplitude
dip=-30;%dip in degrees, positive is down to the right. Try different dips. You might want to change xnot as well
v=2000;
dx=5;
dt=.004;
tmax=1;
x=0:dx:xmax;
t=0:dt:tmax;

filtparms=[5,10,fmax,fmax+20]; %filter specification Ormsby style

%diffraction coordinates
if(dip>0)
    ixnot=round(xnot/dx)+1;
    xd=x(ixnot:ievery:end);
    zd=(xd-xnot)*tand(dip);
    td=2*zd/v;
else
    ixend=round(xnot/dx)+1;
    xd=x(1:ievery:ixend);
    zd=(xd-xnot)*tand(dip);
    td=2*zd/v;
end

%ok install diffractions
seis=zeros(length(t),length(x));
for k=1:length(xd)
    if(between(z1,z2,zd(k)))
        seis=event_hyp(seis,t,x,td(k),xd(k),v,1);
    end
end

%filter
seis=filtf(seis,t,[filtparms(2) filtparms(2)-filtparms(1)],[filtparms(3) filtparms(4)-filtparms(3)]);

plotimage(seis,t,x);
title(['Dipping reflector, dip= ' int2str(dip) ' degrees'])
h=line(xd,td,'color','r');
set(h,'linestyle',':','linewidth',2)
figure
plot(xd,zd);flipy
title('diffraction locations')
grid
axis equal

%fkmig
params=nan*ones(1,13);
params(1)=fmax;
params(3)=dipmax;
[seism,tm,xm]=fkmig(seis,t,x,v,params);
plotimage(seism,tm,xm)
title('Migrated dipping reflector')
h=line(xd,td,'color','r');
set(h,'linestyle',':','linewidth',2)
% spectra
[fk,f,k]=fktran(seis,t,x);
[fkm,fm,km]=fktran(seism,tm,xm);
plotimage(abs(fkm),fm,km)
title('Migrated spectrum')
xlabel('k_x m^{-1}');ylabel('f s^{-1}');
plotimage(abs(fk),f,k)
title('Unmigrated spectrum')
xlabel('k_x m^{-1}');ylabel('f s^{-1}');
%draw evanescent boundaries
fbdy=v*abs(k)/2;
ind=find(fbdy<=max(f));
h=line(k(ind),fbdy(ind),'color','r');
legend(h,'evanescent boundary')
%% #11 horizontal reflector with a local amplitude anomaly
clear
fmax=75;
dipmax=80;
ievery=1;%put a diffraction every this many traces
xmax=2000;
zmax=1000;
dx=5;
zr=zmax/3;%reflector depth;
xwid=10*dx;%width of amplitude anomaly
amp=1.1;%amplitude of the anomaly. Everything else has amplitude 1.0
dz=3*dx;%anomaly will be at depth zr+dz
xend=.1*xmax;%width of no  reflector at edges of section
xnot=xmax/2;%anomaly centered here
v=2000;
dx=5;
dt=.004;
tmax=2*zmax/v;
x=0:dx:xmax;
t=0:dt:tmax;

filtparms=[5,10,fmax,fmax+20]; %filter specification Ormsby style

%diffraction coordinates. Start at xnot and go both directions to ensure
%that we have sampled the anomaly symmetrically
%
inot=round(xnot/dx)+1;
i1=round(xend/dx)+1;
i2=round((xnot-xwid*.5)/dx)+1;
i3=round((xnot+xwid*.5)/dx)+1;
i4=round((xmax-xend)/dx)+1;

xd_left=x(inot:-ievery:i1);
xd_right=x(inot:ievery:i4);
td_left=(2*zr/v)*ones(size(xd_left));
ind=between(xnot,xnot-.5*xwid,xd_left,2);
td_left(ind)=(2*(zr+dz)/v)*ones(size(ind));
td_right=(2*zr/v)*ones(size(xd_right));
ind=between(xnot,xnot+.5*xwid,xd_right,2);
td_right(ind)=(2*(zr+dz)/v)*ones(size(ind));

%ok install diffractions
seis=zeros(length(t),length(x));
for k=1:length(xd_left)
    if(between(xnot,xnot-.5*xwid,xd_left(k),2));
        seis=event_hyp(seis,t,x,td_left(k),xd_left(k),v,amp);
    else
        seis=event_hyp(seis,t,x,td_left(k),xd_left(k),v,1);
    end
end
for k=2:length(xd_right)
    if(between(xnot,xnot+.5*xwid,xd_right(k),2));
        seis=event_hyp(seis,t,x,td_right(k),xd_right(k),v,amp);
    else
        seis=event_hyp(seis,t,x,td_right(k),xd_right(k),v,1);
    end
end

%filter
seis=filtf(seis,t,[filtparms(2) filtparms(2)-filtparms(1)],[filtparms(3) filtparms(4)-filtparms(3)]);


%fkmig
params=nan*ones(1,13);
params(1)=fmax;
params(3)=dipmax;
[seism,tm,xm]=fkmig(seis,t,x,v,params);

data1=seisplotfk(seis,t,x,'Before migration');
data2=seisplotfk(seism,t,x,'After migration');

% plotimage(seism,tm,xm)
% title('Migrated horizontal reflector with anomaly')
% td=2*(zr+dz)/v;
% h=line([xnot-.5*xwid,xnot+.5*xwid],[td td],'color','r');
% set(h,'linestyle','-','linewidth',2)
% legend(h,'anomaly')
% % spectra
% [fk,f,k]=fktran(seis,t,x);
% [fkm,fm,km]=fktran(seism,tm,xm);
% plotimage(abs(fkm),fm,km)
% title('Migrated spectrum')
% xlabel('k_x m^{-1}');ylabel('f s^{-1}');
% plotimage(abs(fk),f,k)
% title('Unmigrated spectrum')
% xlabel('k_x m^{-1}');ylabel('f s^{-1}');
% %draw evanescent boundaries
% fbdy=v*abs(k)/2;
% ind=find(fbdy<=max(f));
% h=line(k(ind),fbdy(ind),'color','r');
% legend(h,'evanescent boundary')
% 
% plotimage(seis,t,x);
% title(['Horizontal reflector, anomaly width ' int2str(xwid) ', amplitude=' num2str(amp) ', depth shift=' num2str(dz)])
% 
% h=line([xnot-.5*xwid,xnot+.5*xwid],[td td],'color','r');
% set(h,'linestyle','-','linewidth',2)
% legend(h,'anomaly')
%% #12
% model a popup reef structure. flanks are always 45 degrees
%
%                           ------------------
%                          .                  .
%                         .                    .
%                        .                      .
%------------------------                        -------------------------
%
%
clear
fmax=75;
dipmax=60;
v=2000;dx=5;dt=.004;%basic model parameters
xmax=3000;
x=0:dx:xmax;%x axis

zb=1000;%depth of base of structure;
zt=800;%depth of top of structure
xw=xmax/3;%width of structure;
xcntr=max(x)/2;%center of structure
amp=1;
tmax=1.5*2*zb/v;

t=0:dt:tmax;%t axis

ievery=1;
seis1=zeros(length(t),length(x));%allocate seismic matrix
seis1=event_diph2(seis1,t,x,v,0,xcntr-.5*xw,zb,ievery,0,amp);
disp('Working ...')
seis1=event_diph2(seis1,t,x,v,xcntr-.5*xw,xcntr-.5*xw+(zb-zt)*tand(45),zb,ievery,-45,amp);
disp('Working ...')
seis1=event_diph2(seis1,t,x,v,xcntr-.5*xw+(zb-zt)*tand(45),xcntr+.5*xw-(zb-zt)*tand(45),zt,ievery,0,amp);
disp('Working ...')
seis1=event_diph2(seis1,t,x,v,xcntr+.5*xw-(zb-zt)*tand(45),xcntr+.5*xw,zt,ievery,45,amp);
disp('Still at it ...')
seis1=event_diph2(seis1,t,x,v,xcntr+.5*xw,max(x),zb,ievery,0,amp);
disp('Wrapping up ...')
[w,tw]=ricker(dt,40,.2);%make ricker wavelet
seis1=sectconv(seis1,t,w,tw);%apply wavelet

%fk migrations
params=nan*ones(1,13);
params(1)=fmax;
params(3)=dipmax;
[seis1mig,tmig,xmig]=fkmig(seis1,t,x,v,params);
plotimage(seis1mig,tmig,xmig)
title(['FK migration of reef with every ' int2str(ievery) 'th diffraction'])

% spectra
[fk,f,k]=fktran(seis1,t,x);
[fkm,fm,km]=fktran(seis1mig,tmig,xmig);
plotimage(abs(fkm),fm,km)
title('Migrated spectrum')
xlabel('k_x m^{-1}');ylabel('f s^{-1}');
plotimage(abs(fk),f,k)
title('Unmigrated spectrum')
xlabel('k_x m^{-1}');ylabel('f s^{-1}');
%draw evanescent boundaries
fbdy=v*abs(k)/2;
ind=find(fbdy<=max(f));
h=line(k(ind),fbdy(ind),'color','r');
legend(h,'evanescent boundary')

plotimage(seis1,t,x);
title(['Popup structure showing every ' int2str(ievery) 'th diffraction hyperbola'])