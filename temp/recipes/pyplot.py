




import matplotlib.pyplot as plt
plt.imshow(time.transpose(),cmap=plt.cm.jet, interpolation='none', extent=[0,5,5,0]);
plt.savefig('test_1_constv.png',format='png',dpi=300,bbox_inches='tight', pad_inches=0)
plt.colorbar(orientation='horizontal',shrink=0.6,label='Traveltime (s)');

fig = plt.figure(figsize=(16, 8))
ax=plt.subplot(5,2,1)
plt.imshow(cmpn.reshape(100,500,order='F'),cmap='jet',clim=(-0.2, 0.2),aspect=0.8)
plt.title('Raw data',color='k');ax.set_xticks([]);ax.set_yticks([]);
ax=plt.subplot(5,2,3)
plt.imshow(dipi.reshape(100,500,order='F'),cmap='jet',clim=(-1,1),aspect=0.8)
plt.title('Iline slope',color='k');ax.set_xticks([]);ax.set_yticks([]);
ax=plt.subplot(5,2,5)
plt.imshow(dipx.reshape(100,500,order='F'),cmap='jet',clim=(-1,1),aspect=0.8)
plt.title('Xline slope',color='k');ax.set_xticks([]);ax.set_yticks([]);
ax=plt.subplot(5,2,7)
plt.imshow(cmpn_d1.reshape(100,500,order='F'),cmap='jet',clim=(-0.2, 0.2),aspect=0.8)
plt.title('Filtered (SOMEAN)',color='k');ax.set_xticks([]);ax.set_yticks([]);
ax=plt.subplot(5,2,9)
plt.imshow((cmpn-cmpn_d1).reshape(100,500,order='F'),cmap='jet',clim=(-0.2, 0.2),aspect=0.8)
plt.title('Noise (SOMEAN)',color='k');ax.set_xticks([]);ax.set_yticks([]);



# ax=plt.subplot(5,2,2)
# plt.imshow(cmpn.reshape(100,500,order='F'),cmap='jet',clim=(-0.2, 0.2),aspect=0.8)
# plt.title('Raw data',color='k');ax.set_xticks([]);ax.set_yticks([]);
# ax=plt.subplot(5,2,4)
# plt.imshow(dipi.reshape(100,500,order='F'),cmap='jet',clim=(-1,1),aspect=0.8)
# plt.title('Iline slope',color='k');ax.set_xticks([]);ax.set_yticks([]);
# ax=plt.subplot(5,2,6)
# plt.imshow(dipx.reshape(100,500,order='F'),cmap='jet',clim=(-1,1),aspect=0.8)
# plt.title('Xline slope',color='k');ax.set_xticks([]);ax.set_yticks([]);
ax=plt.subplot(5,2,8)
plt.imshow(cmpn_d2.reshape(100,500,order='F'),cmap='jet',clim=(-0.2, 0.2),aspect=0.8)
plt.title('Filtered (SOMF)',color='k');ax.set_xticks([]);ax.set_yticks([]);
ax=plt.subplot(5,2,10)
plt.imshow((cmpn-cmpn_d2).reshape(100,500,order='F'),cmap='jet',clim=(-0.2, 0.2),aspect=0.8)
plt.title('Noise (SOMF)',color='k');ax.set_xticks([]);ax.set_yticks([]);
plt.savefig('test_pyseistr_somf3d.png',format='png',dpi=300)

plt.show()

## Important one
plt.gca().set_yticks([]);

plt.figure;
plt.plot(lons,lats,'sr');
plt.plot(lons0,lats0,'pk');
plt.show()

plot(x, y, 'go--', linewidth=2, markersize=12)
plot(x, y, color='green', marker='o', linestyle='dashed',
     linewidth=2, markersize=12)
     
     

e.write("test.qml",format="QUAKEML")


## in plot.waveform

	ntr=len(st)
	
	fig = plt.figure(figsize=(6, 8))
	
	for ii in range(ntr):
		ax = plt.subplot(ntr,1,ii+1)
		nt=len(st[ii].data);
		twin=(nt-1)*1.0/st[ii].stats.sampling_rate;
		staname=st[ii].stats.network+'.'+st[ii].stats.station
		t=np.linspace(0,twin,nt)
		plt.plot(t,st[ii].data,color='k',label = st[ii].id, linewidth = 1, markersize=1)
		
		if ii==0 and titleoff != 1:
			plt.title(st[ii].stats.starttime,fontsize='large', fontweight='normal')
		ax.legend(loc='lower right', fontsize = 10/(ntr/4))
		if ii==ntr-1:
			if axoff == 1:
				plt.setp(ax.get_xticklabels(), visible=False)
			else:
				plt.setp(ax.get_xticklabels(), visible=True)
				ax.set_xlabel("Time (s)",fontsize='large', fontweight='normal')
		else:
			plt.setp(ax.get_xticklabels(), visible=False)
		ax.set_xlim(xmin=0)
		ax.set_xlim(xmax=t[-1])
		ymin, ymax = ax.get_ylim()
		
		if ayoff:
			plt.setp(ax.get_yticklabels(), visible=False)
			
		if picks is not None:
			if picks[ii]['P'] is not None:
				tp=picks[ii]['P']-st[ii].stats.starttime
				plt.vlines(tp, ymin, ymax, color = 'r', linewidth = 1) #for P
			
			if picks[ii]['S'] is not None:
				ts=picks[ii]['S']-st[ii].stats.starttime
				plt.vlines(ts, ymin, ymax, color = 'g', linewidth = 1) #for S
# 		plt.text(2.5,(ymin+ymax)/2,staname,fontsize=12,color='g')

		if ptime is not None:
			tp=ptime[ii]-st[ii].stats.starttime
			plt.vlines(tp, ymin, ymax, color = 'r', linewidth = 1) #for P
			
		if stime is not None:
			ts=stime[ii]-st[ii].stats.starttime
			plt.vlines(ts, ymin, ymax, color = 'g', linewidth = 1) #for P
			
	if figname is not None:
		plt.savefig(figname,**kwargs)
		
		
## Read binary file
fid=open("cchirps.bin","rb");
din = np.fromfile(fid, dtype = np.float32, count = 512).reshape([512,1],order='F')

## Save as binary files
dn=np.float32(dn)
fid = open ("syn3d_dn.bin", "wb") #binary file format, int
fid.write(dn.flatten(order='F'))

d0=np.float32(d0)
fid = open ("syn3d_dc.bin", "wb") #binary file format, int
fid.write(d0.flatten(order='F'))


## compare with matlab
import scipy
from scipy import io
datas = {"d0":d0,"dc":dc,"mask":mask,"dn": dn, "d1": d1, "noi1": noi1, "d2":d2, "noi2":noi2}
scipy.io.savemat("datas3d.mat", datas)

plt.gca().set_ylim(ymin=0,ymax=40);
plt.gca().invert_yaxis();



## Histogram
plt.figure;
plt.hist(deps,10,label='EQCCT',color='b')
plt.hist(deps2,10,label='Catalog',color='g')
plt.gca().set_xlim(xmin=0,xmax=15);
plt.gca().legend(loc='lower right');
plt.gca().set_ylabel("Count",fontsize='large', fontweight='normal')
plt.gca().set_xlabel("Depth (km)",fontsize='large', fontweight='normal')
plt.savefig('continuous_dep_hist.png',format='png',dpi=300)
plt.show() 


## 3D transpose
np.transpose(a, (1, 0, 2)).shape


## specify position of colorbar
#[lower left x, lower left y, upper right x, upper right y] of the desired colorbar:
dat_coord = [-1.5,1.5,-0.5,1.75]
#transform the two points from data coordinates to display coordinates:
tr1 = ax.transData.transform([(dat_coord[0],dat_coord[1]),(dat_coord[2],dat_coord[3])])
#create an inverse transversion from display to figure coordinates:
inv = fig.transFigure.inverted()
tr2 = inv.transform(tr1)
#left, bottom, width, height are obtained like this:
datco = [tr2[0,0], tr2[0,1], tr2[1,0]-tr2[0,0],tr2[1,1]-tr2[0,1]]
#and finally the new colorabar axes at the right position!
cbar_ax = fig.add_axes(datco)
#the rest stays the same:
clevs = [0, 1 , 2]
cb1 = plt.colorbar(hdl, cax=cbar_ax, orientation='horizontal', ticks=clevs)

plt.show()

print(plt.gca().get_xlim());

from matplotlib.ticker import FormatStrFormatter
plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.gca().invert_yaxis()

