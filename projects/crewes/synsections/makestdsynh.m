function [seis,t,x]=makestdsynh(dt,dx,tmax,xmax,v,w,tw)
% MAKESTDSYNH: Make a diffraction synthetic to demo migration codes
%
% [seis,t,x]=makestdsynh(dt,dx,tmax,xmax,v,w,tw)
%
% Make standard migration synthetic. This is a constant velocity 2D synthetic section used to
% demonstrate migration effects. This uses hyperbolic superposition whereas MAKESTDSYN does
% not. So this demo's resolution better. Ths synthetic includes four right-dipping events with dips
% of 20, 40, 60 and 80 degrees, and four left dipping events with the same dips; two noise
% bursts, ten point diffractors, and a horizontal reflector with 5 holes. The two noise spikes
% are at t=.167*tmax,x=xmax/3 and t=.75*tmax,x=2*xmax/3 and have amplitudes 20 times larger
% than the maximum on the horizontal reflector.
%
% dt ... desired time sample rate
% ******* default = .002 *******
% dx ... desired space sample rate
% ******* default =5 ********
% tmax ... maximum time desired in seconds
% ******* default = 2.0 ********
% xmax ... length of section in physical units
% ******* default = 4000 *********
% v ... velocity
% ******* default = 2500 *******
% w ... wavelet
% ******* default 40 Hx Ricker *******
% tw ... wavelet time coordinate
% ******* default is -.1:dt:.1 *******
%
% seis ... synthetic seismic section
% t ... time coordinate vector for seis
% x ... space coordinate vector for seis
% 
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE


if(nargin==6)
    error('if you specify w then you must also give tw');
end
if(nargin<1)
    dt=.002;
end
if(nargin<2)
    dx=5;
end
if(nargin<3)
    tmax=2;
end
if(nargin<4)
    xmax=4000;
end
if(nargin<5)
    v=2500;
end
if(nargin<6)
    [w,tw]=ricker(dt,40,.2);
end

%v=v/2;%exploding reflector synthetic

x1=xmax*4/30;%starting x of right facing dips
x2=xmax*15/30;%ending x of right facing dips
x3=xmax*15/30;%ending x of left facing dips
x4=xmax*26/30;%ending x of left facing dips
%x4=xmax;%ending x of left facing dips

th1=20; th2=40; th3=60; th4=80;

x=0:dx:xmax;
t=0:dt:tmax;

nx = length(x);
nt = length(t);

seis=zeros(nt,nx);
z1=1;%starting depth
%dipping events
disp('eight dipping events to build')
seis= event_diph(seis,t,x,v,x1,x2,z1,th1,1);
disp('completed first event')
seis= event_diph(seis,t,x,v,x1,x2,z1,th2,1);
disp('completed second event')
seis= event_diph(seis,t,x,v,x1,x2,z1,th3,1);
disp('completed third event')
seis= event_diph(seis,t,x,v,x1,x2,z1,th4,1);
disp('completed fourth event')

seis= event_diph(seis,t,x,v,x4,x3,z1,-th1,1);
disp('completed fifth event')
seis= event_diph(seis,t,x,v,x4,x3,z1,-th2,1);
disp('completed sixth event')
seis= event_diph(seis,t,x,v,x4,x3,z1,-th3,1);
disp('completed seventh event')
seis= event_diph(seis,t,x,v,x4,x3,z1,-th4,1);
disp('completed eighth event')

%flat event with holes
t5=tmax/1.5;
z5=t5*v/2;
width=20*dx;
xh01=xmax/20;xh02=xh01+width;
xh1=xmax/4;xh2=xh1+width;
xh3=xmax/2-.5*width;
seis=event_diph(seis,t,x,v,0,xh01,z5,0,1);
seis=event_diph(seis,t,x,v,xh02,xh1,z5,0,1);
seis=event_diph(seis,t,x,v,xh2,xh3,z5,0,1);
xh4=xh3+width;xh5=3*xmax/4-width;
seis=event_diph(seis,t,x,v,xh4,xh5,z5,0,1);
xh6=xh5+width;xh7=xmax-xh01-width;
seis=event_diph(seis,t,x,v,xh6,xh7,z5,0,1);
seis=event_diph(seis,t,x,v,xh7+width,xmax,z5,0,1);


%spikes
ts1=tmax*.167;ts2=tmax*.75;
seis=event_spike(seis,t,x,ts1,xmax/3,100);
seis=event_spike(seis,t,x,ts2,2*xmax/3,100);

%hyperbolas
th1=tmax*.6/1.5;th2=tmax*1.2/1.5;
seis=event_hyp(seis,t,x,th1,xmax/8,v,1);
seis=event_hyp(seis,t,x,th1,xmax/4,v,1);
seis=event_hyp(seis,t,x,th1,xmax/2,v,1);
seis=event_hyp(seis,t,x,th1,3*xmax/4,v,1);
seis=event_hyp(seis,t,x,th1,7*xmax/8,v,1);
seis=event_hyp(seis,t,x,th2,xmax/8,v,1);
seis=event_hyp(seis,t,x,th2,xmax/4,v,1);
seis=event_hyp(seis,t,x,th2,xmax/2,v,1);
seis=event_hyp(seis,t,x,th2,3*xmax/4,v,1);
seis=event_hyp(seis,t,x,th2,7*xmax/8,v,1);


disp('hang on while I apply the wavelet')
seis=sectconv(seis,t,w,tw);

if(nargout==0)
    seisplot(seis,t,x)
    pititle({'The Standard Section with default parameters','made with hyperbolic superposition'})
    seis=[];
    t=[];
    x=[];
end