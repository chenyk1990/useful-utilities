%test decon_workbench
%% #1 stationary test
global MATCHALG
MATCHALG='matchT';
%MATCHALG='matchs';
%make some input data
%s2n=inf;
s2n=2;%signal to noise 
dname='Stationary synthetic data';
dt=.002; %sample rate
tmax=2;%record length
fdom=30;%wavelet dominant frequency
tlen=.5;%wavelet length
ntraces=10;
[r,t]=reflec(tmax,dt,.1,3,pi);
[w,tw]=wavemin(dt,fdom,tlen);
s=convm(r,w);
if(~isinf(s2n))
    s=s+rnoise(s,s2n);
end
seis=s(:,ones(1,ntraces));
% [wr,tw]=ricker(dt,60,tlen);
ind=near(t,.3,1.7);
%sref=convz(r(ind),wr);
sref=r(ind);
tref=t(ind);

x=1:ntraces;
xref=0;
decon_workbench(seis,t,x,sref,tref,xref,'title',[dname ', s2n=' num2str(s2n)])

%% #2 nonstationary test
global MATCHALG
MATCHALG='matchT';
%MATCHALG='matchs';
dname='Nonstationary synthetic data';
%s2n=inf;%signal to noise
s2n=4;
dt=.002; %sample rate
tmax=2;%record length
fdom=20;%wavelet dominant frequency
tlen=.5;%wavelet length
Q=70;%Q value
ntraces=10;
[r,t]=reflec(tmax,dt,.1,3,pi);
[w,tw]=wavemin(dt,fdom,tlen);
qmat=qmatrix(Q,t,w,tw,3);
s=qmat*r;
if(~isinf(s2n))
    s=s+rnoise(s,s2n);
end
seis=s(:,ones(1,ntraces));
% [wr,tw]=ricker(dt,60,tlen);
ind=near(t,.3,1.7);
% sref=convz(r(ind),wr);
sref=r(ind);
tref=t(ind);
x=1:ntraces;
xref=0;
decon_workbench(seis,t,x,sref,tref,xref,'title',[dname ', s2n=' num2str(s2n)])

%% #3 test on a panel of real data
global MATCHALG
MATCHALG='matchT';
%MATCHALG='matchs';
dname='Real data';
%load real data test panel
load decon_workbench_real_test_panel

decon_workbench(spanel(:,1:2:end),t,xpanel(1:2:end),sref(ind),tref(ind),618,'title',dname);

