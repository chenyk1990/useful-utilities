%% Given a seismic section with a 'single' live sample show the corresponding migration wavefront 
%
% Change the x and t coordinates of the live point and observe the change in the corresponding
% wavefronts. 
%

xmax=2000;%line length
v=2000;
tmax=1;
xnot=xmax/2;%x coordinate of live point
tnot=tmax/2;%t coordinate of the live point

dx=5;
dt=.002;
x=0:dx:xmax;
t=0:dt:tmax;

tx=sqrt(tnot.^2+4*((x-xnot)/v).^2);
filtparms=[5,10,80,100]; %filter specification Ormsby style
seis=zeros(length(t),length(x));

seism=event_wavefront(seis,t,x,tnot,xnot,v,1);

it=near(t,tnot);
ix=near(x,xnot);
seis(it(1),ix(1)-3:ix(1)+3)=1;

zmig=v*t/2;
%filter
seis=filtf(seis,t,[filtparms(2) filtparms(2)-filtparms(1)],[filtparms(3) filtparms(4)-filtparms(3)]);
seism=filtf(seism,t,[filtparms(2) filtparms(2)-filtparms(1)],[filtparms(3) filtparms(4)-filtparms(3)]);

datar=seisplottwo(seis,t,x,'Seismic section with one live sample',seism,zmig,x,...
    {'The corresponding wavefront circle','aka the depth picture'});
axes(datar{2})
axis equal
ylim([0 zmig(end)])

%% place wavefronts along a hyperbolic path
%
% This is a geometric demonstration of focussing a hyperbolic diffraction by wavefront
% replacment. Essentiall each point in an unmigrated image is replaced by a circular wavefront
% (2D) or spherical wavefront (3D), where the amplitude along the wavefront is the amplitude of
% the point being replace. Thus each point is speread along a wavefront. If tp is the time of
% the sample being replaced, then the radius of the wavefront is r=tp*v/2. If the seismic
% section contains only a single diffraction curve, then only those wavefronts from points on
% the diffraction are non-zero. The result is then shown here. Run this for ievery = 20 and
% then ievery=1. View the results at different clip levels and convince yourself that this
% process of wavefront replacement focusses the energy on a diffraction curve.
%
ievery=1;%place a wavefront every this many traces
xmax=2000;
tmax=1;
xnot=1000;%x coordinate of apex of hyperbola
tnot=.25;%apex time of the hyperbola


v=2000;
dx=5;
dt=.002;
x=0:dx:xmax;
t=0:dt:tmax;

tx=sqrt(tnot.^2+4*((x-xnot)/v).^2);
filtparms=[5,10,80,100]; %filter specification Ormsby style
seis=zeros(length(t),length(x));

for k=1:ievery:length(x)
    seis=event_wavefront(seis,t,x,tx(k),x(k),v,1);
end

%filter
seis=filtf(seis,t,[filtparms(2) filtparms(2)-filtparms(1)],[filtparms(3) filtparms(4)-filtparms(3)]);

seisplot(seis,t,x);
h=line(x,tx,'color','r','linestyle','none','marker','.');
pititle('Wavefront circles focus a hyperbolic diffraction')

%% migrate a diffraction by wavefont superposition

ieveryt=10;%(migrating)place a wavefront every this many time samples
ieveryx=1;%(migrating)place a wavefront every this many traces
xmax=2000;%line length
zmax=1000;%maximum depth
xnot=xmax/2;%x coordinate of the diffractor
znot=zmax/2;%z coordinate of the diffractor
tnot=2*znot/v;%apec of the hyperbola
v=2000;%velocity
dx=5;%grid size
dt=.004;%time sample rate
tmax=zmax*2/v;%maximum time

x=0:dx:xmax;
t=0:dt:tmax;
z=0:dx:zmax;
filtparms=[5,10,80,100]; %filter specification Ormsby style
%Create the diffraction
seis=zeros(length(t),length(x));
seis=event_hyp(seis,t,x,tnot,xnot,v,1);
seis=filtf(seis,t,[filtparms(2) filtparms(2)-filtparms(1)],[filtparms(3) filtparms(4)-filtparms(3)]);
seismig=wavefrontmig(seis,t,x,v,xmax,filtparms,[ieveryt ieveryx]);
zmig=t*v/2;

datar=seisplottwo(seis,t,x,{'A single diffraction',['velocity = ' num2str(v)]},...
    seismig,zmig,x,{'Wavefront migration of the diffraction',['ieveryt=' int2str(ieveryt)...
    ', ieveryx=' int2str(ieveryx)]});
axes(datar{2})
axis equal
ylim([0 zmig(end)])
%% illustrate need for balanced amplitides migrate a diffraction by wavefont superposition

ieveryt=10;%(migrating)place a wavefront every this many time samples
ieveryx=1;%(migrating)place a wavefront every this many traces
pctamp=200;%pecent of amplitude disturbance
xmax=2000;%line length
zmax=1000;%maximum depth
xnot=xmax/2;%x coordinate of the diffractor
znot=zmax/2;%z coordinate of the diffractor
tnot=2*znot/v;%apec of the hyperbola
v=2000;%velocity
dx=5;%grid size
dt=.004;%time sample rate
tmax=zmax*2/v;%maximum time

x=0:dx:xmax;
t=0:dt:tmax;
z=0:dx:zmax;
filtparms=[5,10,80,100]; %filter specification Ormsby style
%Create the diffraction
seis=zeros(length(t),length(x));
seis=event_hyp(seis,t,x,tnot,xnot,v,1);
seis=filtf(seis,t,[filtparms(2) filtparms(2)-filtparms(1)],[filtparms(3) filtparms(4)-filtparms(3)]);
%distort trace amplitudes in a random way
arand=.01*pctamp*randn(size(x))+1;
seis2=seis.*arand(ones(size(t)),:);
seismig=wavefrontmig(seis,t,x,v,xmax,filtparms,[ieveryt ieveryx]);
seismig2=wavefrontmig(seis2,t,x,v,xmax,filtparms,[ieveryt ieveryx]);
zmig=t*v/2;

datar=seisplottwo(seismig,zmig,x,{'Migrated diffraction',['velocity = ' num2str(v)]},...
    seismig2,zmig,x,{'Migrated diffration with disturbed amplitudes',['ieveryt=' int2str(ieveryt)...
    ', ieveryx=' int2str(ieveryx)]});
% axes(datar{2})
% axis equal
% ylim([0 zmig(end)])
%% migrate an entire diffraction chart by wavefront superposition
% construct the depth and time pictures of a vertical array of scatterpoints
% What do the asymptotes tell you? What is the macimum possible time dip?

%constant velocity
znot=[40 60 100 150 200 300 400 600 800 1000];%these are the scatterpoint depths
xmax=2000;%line length
zmax=1000;%max depth
xnot=xmax/2;%x coordinate of the diffractor

v=2000;
dx=2;
dt=.002;
tmax=zmax*2/v;
x=0:dx:xmax;
t=0:dt:tmax;
z=0:dx:zmax;
filtparms=[5,10,80,100]; %filter specification Ormsby style
%Create the diffraction

seis=zeros(length(t),length(x));
for k=1:length(znot)
    tnot=2*znot(k)/v;
    seis=event_hyp(seis,t,x,tnot,xnot,v,1);
end
seis=filtf(seis,t,[filtparms(2) filtparms(2)-filtparms(1)],[filtparms(3) filtparms(4)-filtparms(3)]);
% create the depth section
depth=zeros(length(z),length(x));
for k=1:length(znot)
    ix=near(x,xnot);
    iz=near(z,znot(k));
    depth(iz,ix)=1;
end
g=gaus_radial(dx,5*dx);
depth=conv2(depth,g,'same');

datar=seisplottwo(depth,z,x,'Array of scatterpoints in depth',seis,t,x,{'A diffraction chart',['velocity = ' num2str(v)]});
axes(datar{1})
axis equal
ylim([0 zmax])
axes(datar{2})
ta=2*abs(x-xnot)/v;
h=line(x,ta,'linestyle','--','color','r');
ylabel('two-way time (seconds)')
legend(h,'hyperbolic asymptotes')
%% migrate a dipping reflector by wavefront superposition
%
% Here we run the program wavefrontmig which migrates by wavefront replacement. This is
% actually a form of Kirchhoff migration. There are 3 ievery parameters now, ievery is used in
% modelling and you are familiar with it. ieveryt and ieveryx are used in migration and
% determine the spacing of wavefronts in t and x. If ieveryt=10 and ieveryx=1 then every tenth
% sample of every trace is replaced by a wavefrot.
%
% Try this and notice how the unmigrated image is formed as a tangent to the diffraction curves
% while the migrated image is a tangent to the wavefronts.
%
ievery=1;%(modelling) place a diffractor every this many traces. A value of 1 give the full picture
ieveryt=1;%(migrating)place a wavefront every this many time samples
ieveryx=1;%(migrating)place a wavefront every this many traces
dip=30;%dip in degrees, positive is down to the right. Try different dips. You might want to change xnot as well
xmax=2000;
zmax=1000;
xnot=500;%x coordinate where the reflector outcrops
z1=100;%first depth with nonzero amplitude
z2=1000;%last depth with nonzero amplitude
v=2000;
dx=10;
dt=.004;
tmax=1;
x=0:dx:xmax;
t=0:dt:tmax;

filtparms=[5,10,60,80]; %filter specification Ormsby style

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
seismig=wavefrontmig(seis,t,x,v,xmax,filtparms,[ieveryt ieveryx]);
zmig=t*v/2;

if(ievery==1)
    titl1={['Unmigrated Dipping reflector, dip= ' int2str(dip) ' degrees'],'Every point replaced'};
else
    titl1={['Unmigrated Dipping reflector, dip= ' int2str(dip) ' degrees'],['Every ' int2str(ievery) 'th point replaced']};
end

datar=seisplottwo(seis,t,x,titl1,seismig,zmig,x,{'Dipping reflector migrated by wavefront superposition',...
    ['ieveryt=' int2str(ieveryt) ', ieveryx=' int2str(ieveryx)]});

axes(datar{1})
h=line(xd,td,'color','r');
set(h,'linestyle','none','linewidth',2,'marker','.');
legend(h,'scatterpoints')
axes(datar{2})
h=line(xd-2*dx,zd,'color','r');
set(h,'linestyle','none','linewidth',2,'marker','.');
legend(h,['scatterpoints shifted ' num2str(2*dx) 'm to the left'])
axis equal
ylim([0 zmig(end)])

%% migrate a dip fan by wavefront superposition

ievery=1;%(modelling) place a diffractor every this many traces. A value of 1 give the full picture
ieveryt=1;%(migrating)place a wavefront every this many time samples
ieveryx=1;%(migrating)place a wavefront every this many traces
xmax=2000;
tmax=1;
xnot=500;%x coordinate of outcrop (t=0) of linear event
dips=0:15:60;%dips in degrees of the events
tnot=.2;


v=2000;
dx=10;
dt=.004;
x=0:dx:xmax;
t=0:dt:tmax;
filtparms=[5,10,60,80]; %filter specification Ormsby style
seis=zeros(length(t),length(x));
diplbl='dips: ';
for k=1:length(dips)
    x1=xnot;
    x2=.75*xmax;
    seis=event_diph(seis,t,x,v,x1,x2,(tnot+(k-1)*.1*tnot)*v/2,dips(k),[1 1]);
    if(k~=length(dips))
        diplbl=[diplbl ' ' num2str(dips(k)) ','];
    else
        diplbl=[diplbl ' ' num2str(dips(k)) ' degrees'];
    end
end

%filter
seis=filtf(seis,t,[filtparms(2) filtparms(2)-filtparms(1)],[filtparms(3) filtparms(4)-filtparms(3)]);
seismig=wavefrontmig(seis,t,x,v,xmax,filtparms,[ieveryt ieveryx]);
zmig=t*v/2;

datar=seisplottwo(seis,t,x,'Unmigrated dip fan',seismig,zmig,x,{'migrated dip fan',diplbl});
axes(datar{2})
axis equal
ylim([0 zmax])