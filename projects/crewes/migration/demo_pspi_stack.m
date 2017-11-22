%% Make the input data
% this takes 20 minutes. But the result is saved to disk so you only need to do it once.

%exploding reflector model of thrust

modelname='thrust model';
dx=2.5;
[velt,xt,zt]=thrustmodel(dx);

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

%%
%pspi depth migration
%this if block looks for the saved dataset from the previous cell
if(~exist('seisex','var'))
    if(exist('thrust_explode.mat','file'))
        load thrust_explode
    else
        disp('Please run the first cell of this script')
    end
end

%establish extrapolation checkpoints (see help for pspi_stack)
zcheck=0:50:zt(end);
dt=tt(2)-tt(1);
dx=xt(2)-xt(1);
%velocity model is sampled every 2.5 m, the next line allows a 5m depth step.
%ideally, we showld average evergy two depths into one (i.e. resample the
%model) but that is a detail.
zmig=zt(1:2:end);
[seismig,checkinfo]=pspi_stack(seisex,tt,xt,velt(1:2:end,:),xt,zmig,[0 inf],zcheck);
plot_pspi_checkinfo(checkinfo,'Thrust model')

seisplottwo(seisex,tt,xt,'Thrust model: Unmigrated ERM section',seismig,zmig,xt,'PSPI depth migration');

imagesc(xt,zt,velt);colorbar
title('Thrust velocity model')

%%
%kirchhoff time migration
if(~exist('seisex','var'))
    if(exist('thrust_explode.mat','file'))
        load thrust_explode
    else
        disp('Please run the first cell of this script')
    end
end
%make an rms velcity model
tmax=max(tt);
dt=tt(2)-tt(1);
[velrms,tv]=vzmod2vrmsmod(velt,zt,dt,tmax,2);

figure
imagesc(xt,tv,velrms);colorbar
title('RMS velocity model')
[seistmig,tmig,xmig]=kirk(seisex,velrms,tt,xt);

plotimage(seistmig,tmig,xmig);
title('Kirkhhoff time migration')


%%
% Demo time and depth migration on a flat section.
% Here we make a section with flat events and migrate it with the thrust
% velocity model. The results are "pleasing" with time migration but not
% with depth migration. This demonstrates the built-in bias of time
% migration to preserve flat (horizontal) events. In actuality, a migration
% such as this is a serious mistake and depth migration, by getting a bad
% result, lets us know that we are wrong. Time migration give the false
% impression that all is well.
%
% Be sure you've run the previous sections first
%
if(~exist('seisex','var'))
    if(exist('thrust_explode.mat','file'))
        load thrust_explode
    else
        disp('Please run the first cell of this script')
    end
end
% make a flat section
dt=t(2)-t(1);
r=reflec(max(tt),dt);
[w,tw]=ricker(dt,40,.2);
s=convz(r,w);
seisflat=s*ones(size(xt));
plotimage(seisflat,tt,xt)
title('flat section')
%time migration 
[flat_tmig,tmig,xmig]=kirk(seisflat,velrms,tt,xt);

plotimage(flat_tmig,tmig,xmig);
title('Kirkhhoff time migration of flat section')

%depth migration
zmig=zt(1:2:end);
flat_dmig=pspi_stack(seisflat,tt,xt,velt(1:2:end,:),xt,zmig);
plotimage(flat_dmig,zmig,x)
title('PSPI depth migration of flat section')