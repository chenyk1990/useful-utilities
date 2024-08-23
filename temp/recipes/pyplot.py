import warnings
warnings.simplefilter("ignore")

import sys
sys.path.append('.')


plt.imshow(cmpn.reshape(100,500,order='F'),cmap='jet',clim=(-0.2, 0.2),aspect='auto')

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

#imshow with axis scale  ax=plt.imshow(self.vp[:,:,ii],extent=[self.minx,self.maxx,self.maxz,self.minz],aspect='auto');
        

datan=np.concatenate(datan,axis=2)

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
     linewidth=2, markersize=12, fillstyle='none') #fillstyle='full'
     
     

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


## text box
def addtext(ax, props):
    ax.text(0.5, 0.5, 'text 0', props, rotation=0)
    ax.text(1.5, 0.5, 'text 45', props, rotation=45)
    ax.text(2.5, 0.5, 'text 135', props, rotation=135)
    ax.text(3.5, 0.5, 'text 225', props, rotation=225)
    ax.text(4.5, 0.5, 'text -45', props, rotation=-45)
    for x in range(0, 5):
        ax.scatter(x + 0.5, 0.5, color='r', alpha=0.5)
    ax.set_yticks([0, .5, 1])
    ax.set_xticks(np.arange(0, 5.1, 0.5))
    ax.set_xlim(0, 5)
    ax.grid(True)


# the text bounding box
bbox = {'fc': '0.8', 'pad': 0}

fig, axs = plt.subplots(2, 1, sharex=True)

addtext(axs[0], {'ha': 'center', 'va': 'center', 'bbox': bbox})


## compare with matlab
import scipy
from scipy import io
datas = {"d0":d0,"dc":dc,"mask":mask,"dn": dn, "d1": d1, "noi1": noi1, "d2":d2, "noi2":noi2}
scipy.io.savemat("datas3d.mat", datas)

datas=scipy.io.loadmat("event_6~16.mat") (dict type)

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

plt.axis('off')

cax = fig.add_axes([0.15,0.9,0.2,0.2])
plt.text(0,0, "a)", fontsize=28, color='k') 
plt.gca().text(-0.15,1,'(a)',transform=plt.gca().transAxes,size=20,weight='normal')
plt.axis('off')
## difference catalog (what is missing)
eids2=[ii for ii in eids if ii not in eids0]#there are 703 events in eids,


if os.path.isdir('./newevents') == False:  
	os.makedirs('./newevents',exist_ok=True)
os.path.isfile('test.dat')


import obspy.core.utcdatetime as utc
utc.UTCDateTime(year, month, day, 00, 00, 00, 000000)+86400)
	

TypeError("Only integers are allowed")

Exception("Sorry, no numbers below zero")


from pylib.texnet import read_events
e=read_events(['texnet2020galz'])

from pylib.texnet import eventlist
e=eventlist(['texnet2020galz'])

from pylib.texnet import cityloc
cityloc('Midland')
	

## add a scalebar
pip install matplotlib-scalebar



ot='2016-10-26T03:15:36.000000Z'
str(ot).replace("-","").replace(":","") 


	from pyseistr import cseis
	import numpy as np
	from matplotlib import pyplot as plt
	plt.imshow(np.random.randn(100,100),cmap=cseis())
	plt.show()
	
		try:
			e = cl.get_events(eventid=eids[ie], includearrivals=True)
			s=get_streams_frompicks(e[0],tbefore=180,tafter=-120);
			s.write(filename='./allnoise/'+eids[ie]+'_1.mseed', format="MSEED")

		except:
			print("Event ID %s is wrong !!!!!!!!"%(eids[ie]))
		else:
			pass
			
#hide legend
plt.hist(mags,40,range=(np.min(mags),np.max(mags)),label=None,color='lightgray',edgecolor='black',log=True)
plt.gca().legend().set_visible(False)

Bbox=plt.gca().get_position()

xmin,ymin=plt.gca().get_position().get_points()[0]
ymax,ymax=plt.gca().get_position().get_points()[1]
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
plt.gca().axis('off')

plt.gcf().add_axes([xmin,ymin,0.2,0.2])


## 3D scatter
import numpy
from numpy.random import rand

xs = rand(10)
ys = rand(10)
zs = rand(10)
amps = rand(10)

import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(xs, ys, zs, c=amps)
plt.show()

%% add patch
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
#define Matplotlib figure and axis
fig, ax = plt.subplots()
#create simple line plot
ax.plot([0, 10],[0, 10])
#add rectangle to plot
ax.add_patch(Rectangle((1, 1), 2, 6, alpha=0.5)) #xy,width,height
#display plot
plt.show()


#save pandas frame to csv
df.to_csv('test.csv')

import pandas as pd
df2=pd.read_csv("test.csv")


## plot loss
print(history.history.keys())
#  "Accuracy"
plt.plot(history.history['acc'])
plt.plot(history.history['val_acc'])
plt.title('model accuracy')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train', 'validation'], loc='upper left')
plt.show()
# "Loss"
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'validation'], loc='upper left')
plt.show()


## pandas
df[df['Week'] == 'Week7'].sort_values(by=['Time']).iloc[0]

df[df['Week'] == 'Week7'].sort_values(by=['Time']).iloc[0:3,3:4]

df[df['Week'] == 'Week7'].sort_values(by=['Time']).iloc[0:3,3:4]['mag']


tmpframe[tmpframe['Mag']==tmpframe['Mag'].max()]


#numpy 
np.array_equal(dataall,dd)

datan=np.load('data.npy')
np.save('data.npy',data)

np.expand_dims(data,1)

#numpy 
#pandas,dataframe
df.reset_index(drop=True)

