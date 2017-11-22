%% make the input data
modelname='thrust model';
dx=2.5;
[velt,xt,zt]=thrustmodel(dx);
zmax=max(zt);
xmax=max(xt);
vlow=min(velt(:));
vmax=max(velt(:));

figure
imagesc(xt,zt,velt);colorbar
colormap(copperud)
title(modelname)
xlabel('distance (m)');
ylabel('depth (m)');
prepfig

dt=.004; %temporal sample rate
dtstep=.0005;
tmax=2*zmax/vlow; %maximum time
fdom=35;
[w,tw]=wavemin(dt,fdom,.5);
[seisex,seisexnw,tt]=afd_explode(dx,dtstep,dt,tmax, ...
 		velt,xt,zeros(size(xt)),w,tw,2);
    
fnameout='thrust_explode';
save(fnameout,'seisex','tt','xt','zt','velt','w','tw')

seisplot(seisex,tt,xt,'Thrust exploding reflector model')

    
%% now an RTM

if(exist('thrust_explode.mat','file'))
    load thrust_explode
    disp('thrust_explode loaded from disk')
end

dtstep=.0005;

tout=0:.05:2;

[seism,xm,zm,snapshots,xms,zms]=rtm_stack(seisex,tt,xt,dtstep,velt,tout);

seisplot(seism,zm,xm,'Post-stack RTM migration of Thrust model')

tout=fliplr(tout);
titles=cell(size(snapshots));
for k=1:length(titles)
    titles{k}=['snapshot at time tout=' num2str(tout(k))];
end

tm=2*zm/3000;
%seismf=butterband(seism,tm,5,0,4,0);
%highpass filter
seismf=filter_stack(seism,tm,5,0,'method','filtf');

seisplottwo(seisex,tt,xt,'Unmigrated',seism,zm,xm,'Migrated');
seisplottwo(seismf,zm,xm,'Migrated filtered',seism,zm,xm,'Migrated');

plotsnaps(snapshots(1:41),velt,xms,zms,'Distance (m)','Depth (m)',titles(1:41),'Thrust model','copperud','Post-Stack RTM');

fnameout='thrust_explode_rtm';
save(fnameout,'seism','xm','zm','snapshots','xms','zms','seisex','tt','xt','velt','zt','tout','titles');