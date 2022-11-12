




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
		
		


