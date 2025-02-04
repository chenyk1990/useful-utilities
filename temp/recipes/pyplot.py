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


plt.colorbar(orientation='horizontal',cax=fig.add_axes([0.37,0.07,0.3,0.02]),shrink=1,label='Traveltime (s)');

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

plt.colorbar(orientation='horizontal',shrink=0.6,label='Traveltime (s)');

plt.show()

print(plt.gca().get_xlim());

from matplotlib.ticker import FormatStrFormatter
plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.gca().invert_yaxis()

plt.axis('off')

string='hanning{0}.pdf'.format(2.0,'3',5)
string='hanning{1}.pdf'.format(2.0,'3',5)

txt = "For only {price:.2f} dollars!"
print(txt.format(price = 49))

plt.title('Clean:Noisy:Denoised:Noise {SNR:.2f}'.format(SNR=snr(data,denoised)))


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

#matplotlib

#remove margin
	import matplotlib.pyplot as plt	
	plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0); 
	plt.imshow(dn,clim=(-1, 1),aspect='auto');plt.axis('off');plt.margins(0, 0);plt.savefig('dn.png')
	plt.imshow(d1,clim=(-1, 1),aspect='auto');plt.margins(0, 0);plt.savefig('d1.png')


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

## 
np.random.seed(1337)
numpy.random.seed(1337)
numpy.random.shuffle()
HFef=asciiread(os.getenv('HOME')+'/chenyk.hfswd/evetmp/eveHFef.dat')
HFef=np.array(HFef)
np.random.seed(20232425)
print(HFef[0:2])
np.random.shuffle(HFef)
print(HFef[0:2])

ind=int(len(HFef)*0.9)
HFef1=HFef[0:ind]
HFef2=HFef[ind:]

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


## deep learning
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Convolution2D, MaxPooling2D
from keras.utils import to_categorical

batch_size = 128
nb_classes = 10
nb_epoch = 12

# the data, shuffled and split between tran and test sets
(X_train, y_train), (X_test, y_test) = mnist.load_data()

# X_train = X_train.reshape(X_train.shape[0], 1, 28, 28) #torch NCHW
# X_test = X_test.reshape(X_test.shape[0], 1, 28, 28)	 #torch NCHW
X_train = X_train.reshape(X_train.shape[0], 28, 28, 1)	 #tf NHWC
X_test = X_test.reshape(X_test.shape[0], 28, 28, 1)		 #tf NHWC

X_train = X_train.astype("float32")
X_test = X_test.astype("float32")
X_train /= 255
X_test /= 255
print('X_train shape:', X_train.shape)
print(X_train.shape[0], 'train samples')
print(X_test.shape[0], 'test samples')

# convert class vectors to binary class matrices
Y_train = to_categorical(y_train, nb_classes)
Y_test = to_categorical(y_test, nb_classes)

model = Sequential()

model.add(Convolution2D(32, 3))
model.add(Activation('relu'))
# model.add(Convolution2D(32, 3))
# model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128))
model.add(Activation('relu'))
model.add(Dropout(0.5))
model.add(Dense(nb_classes))
model.add(Activation('softmax'))

model.compile(loss='categorical_crossentropy', optimizer='adadelta', metrics=['accuracy'])

history=model.fit(X_train, Y_train, batch_size=batch_size, epochs=nb_epoch, verbose=1, validation_data=(X_test, Y_test))
score = model.evaluate(X_test, Y_test, verbose=1)
print('Test score:', score[0])
print('Test accuracy:', score[1])

print(history.history.keys())
import matplotlib.pyplot as plt
#  "Accuracy"
plt.plot(history.history['accuracy'])#sometimes it's history.history['acc']
plt.plot(history.history['val_accuracy']) #sometimes it's history.history['val_acc']?
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

ax.xaxis.set_tick_params(labelsize=20)
ax.yaxis.set_tick_params(labelsize=20)
  

## pandas
df[df['Week'] == 'Week7'].sort_values(by=['Time']).iloc[0]

df[df['Week'] == 'Week7'].sort_values(by=['Time']).iloc[0:3,3:4]

df[df['Week'] == 'Week7'].sort_values(by=['Time']).iloc[0:3,3:4]['mag']


tmpframe[tmpframe['Mag']==tmpframe['Mag'].max()]



#in matlab
a(find(a>0.6))=0 the same as: a(a>0.6)=0
  
#numpy 
np.ma.masked_array (fill some positions with values)
example
a=np.ma.masked_array(
   data=[[5, 8, 17],
         [8, 16, 24],
         [17, 24, 61]],
   mask=[[False, False, False],
         [False, False, False],
         [False, False, True]],fill_value=999999)
    
rng.choice(5, 3, replace=False) #This is equivalent to rng.permutation(np.arange(5))[:3]

np.isfinite(sp_amp) -> MATLAB: isfinite(sp_amp)

np.atleast_2d(3.0)
array([[3.]])

rng.integers  #    
#randomized_vm_ind=rng.integers(low=0, high=num_velocity_models,size=perturbed_origin_depth_km.shape[-1]);
#rand('state',123);randomized_vm_ind=round(rand(size(perturbed_origin_depth_km,2),1)*num_velocity_models);

rng=np.random.default_rng(123) #seed is 123
rng.normal(size=(nmc,num_unique_events)) #Draw random samples from a normal (Gaussian) distribution.

.astype(int)

tmp=(x[idx,np.arange(n2)]==False) #This is very tricky, [x[idx,ii] for ii in range(n2)]
	
ytop1=utop-ptab[:, np.newaxis] 
#dimension: [1xn2] - [n1x1] -> [n1xn2] = Matlab: ones(n1,1)*utop(:)' - ptab(:)*ones(1,n2)
ytop2=utop+ptab[:, np.newaxis] #dimension: [1xn2] - [n1x1] -> [n1xn2]

np.cross(a,b) = cross(a,b)

np.vstack #(Stack arrays in sequence vertically (row wise).)
np.vstack((vmodel_depthvp,vmodel_depthvp[-1,:]) -> Matlab [vmodel_depthvp;vmodel_depthvp(end,:)]

np.hstack #(Stack arrays in sequence horizontally (column wise).)

np.insert #Insert values along the given axis before the given indices.
#np.insert(z,i,z(i)) -> [z(1:i),z(i+1),z(i+1:end)];

np.nansum #Return the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.

if id1.ndim==1:
	id1[np.newaxis,:].shape #create the first axis [1xdim]
	id1=np.repeat(id1[np.newaxis,:],size,axis=0) #spray vector id1/id2 along first axis
t=t[:,:,np.newaxis] = np.expand_dims(t,2)

np.arctan2 = atan2
 
x.argmax -> [~,t]=max(x)
x.argmin -> [~,t]=min(x)

np.diff = diff (matlab)
np.cumsum=cumsum (matlab)

np.isin (ismember(1,a)?)

np.repeat([1,2,3],3,axis=0)
numpy.any([1,2,0,0]) #if any member of the list is >0 

np.arange(0, 5.1, 0.5)
np.arange(10)

np.array_equal(dataall,dd)

datan=np.load('data.npy')
np.save('data.npy',data)

np.expand_dims(data,1)

np.linalg.svd
np.linalg.norm

# a is an numpy array
print(np.isrealobj(a)) # -> True 
print(np.iscomplexobj(a)) # -> False

scalogram2 = np.empty(dout.shape[0:2], dtype=np.complex64)
scalogram2.real = dout[:,:,0]
scalogram2.imag = dout[:,:,1]

np.ma.divide -> Divide arguments element-wise. (Equivalent to ``x1`` / ``x2`` in terms of array-broadcasting.)

#numpy 
#pandas,dataframe
df.reset_index(drop=True)

a = [1,2,3,4,1,2,3,4,1,2,3,4]
a[::2] #every two
a[::3] #every three

#numpy thresholding
tmp=dout[:,:,0]*dout[:,:,0]+dout[:,:,1]*dout[:,:,1]
tmp=tmp/tmp.max()
tmp[tmp<0.05]=0
tmp[tmp>=0.05]=1

#numpy index
dtf0[:,np.linspace(100,300,nf,dtype='int')].reshape([1,n1,nf,1])

#
np.bool -> np.bool_ #module 'numpy' has no attribute 'bool'.

import matplotlib.pyplot as plt
from pyseistr import framebox

x1=200
x2=400
y1=1500
y2=5000

plt.figure(figsize=(25, 16))
plt.subplot(231)
plt.imshow(data[:,:],aspect='auto',vmin=-20,vmax=20,cmap="seismic");
plt.xlabel('Channel',size=16,weight='bold');plt.ylabel('Sample',size=16,weight='bold');
plt.title('Original (115.2 Mb)',size=16, weight='bold');
framebox(x1,x2,y1,y2);
plt.gca().text(-0.17,1,'a)',transform=plt.gca().transAxes,size=20,weight='bold')
plt.gca().xaxis.set_tick_params(labelsize=16)
plt.gca().yaxis.set_tick_params(labelsize=16)
plt.rc('font', weight='bold');

plt.subplot(232)
plt.imshow(outB[:,:],aspect='auto',vmin=-20,vmax=20,cmap="seismic")
plt.xlabel('Channel',size=16,weight='bold');
plt.title('Reconstructed (19.14 Mb)',size=16, weight='bold');
framebox(x1,x2,y1,y2);
plt.gca().text(-0.17,1,'b)',transform=plt.gca().transAxes,size=20,weight='bold')
plt.gca().xaxis.set_tick_params(labelsize=16)
plt.gca().yaxis.set_tick_params(labelsize=16)
plt.rc('font', weight='bold');

plt.subplot(233)
plt.imshow(data[:,:]-outB[:,:],aspect='auto',vmin=-20,vmax=20,cmap="seismic")
plt.xlabel('Channel',size=16,weight='bold');
plt.title('Error',size=16, weight='bold');
framebox(x1,x2,y1,y2);
plt.gca().text(-0.17,1,'c)',transform=plt.gca().transAxes,size=20,weight='bold')
plt.gca().xaxis.set_tick_params(labelsize=16)
plt.gca().yaxis.set_tick_params(labelsize=16)
plt.rc('font', weight='bold');

plt.subplot(234)
plt.imshow(data[y1:y2,x1:x2],aspect='auto',vmin=-100,vmax=100,cmap="seismic",extent=[x1,x2,y2,y1]);
plt.xlabel('Channel',size=16,weight='bold');plt.ylabel('Sample',size=16,weight='bold');#plt.title('Original (115.2 Mb)')
plt.gca().spines['bottom'].set_color('red')
plt.gca().spines['top'].set_color('red') 
plt.gca().spines['right'].set_color('red')
plt.gca().spines['left'].set_color('red')
plt.gca().text(-0.17,1,'d)',transform=plt.gca().transAxes,size=20,weight='bold')
plt.gca().xaxis.set_tick_params(labelsize=16)
plt.gca().yaxis.set_tick_params(labelsize=16)
plt.rc('font', weight='bold');

plt.subplot(235)
plt.imshow(outB[y1:y2,x1:x2],aspect='auto',vmin=-100,vmax=100,cmap="seismic")
plt.xlabel('Channel',size=16,weight='bold');#plt.title('Reconstructed (19.14 Mb)');
plt.gca().spines['bottom'].set_color('red')
plt.gca().spines['top'].set_color('red') 
plt.gca().spines['right'].set_color('red')
plt.gca().spines['left'].set_color('red')
plt.gca().text(-0.17,1,'e)',transform=plt.gca().transAxes,size=20,weight='bold')
plt.gca().xaxis.set_tick_params(labelsize=16)
plt.gca().yaxis.set_tick_params(labelsize=16)
plt.rc('font', weight='bold');

plt.subplot(236)
plt.imshow(data[y1:y2,x1:x2]-outB[y1:y2,x1:x2],aspect='auto',vmin=-100,vmax=100,cmap="seismic")
plt.xlabel('Channel',size=16,weight='bold');#plt.title('Error')
plt.gca().text(-0.17,1,'f)',transform=plt.gca().transAxes,size=20,weight='bold')
plt.gca().spines['bottom'].set_color('red')
plt.gca().spines['top'].set_color('red') 
plt.gca().spines['right'].set_color('red')
plt.gca().spines['left'].set_color('red')
plt.gca().xaxis.set_tick_params(labelsize=16)
plt.gca().yaxis.set_tick_params(labelsize=16)
plt.rc('font', weight='bold');

plt.savefig('forge1.png',dpi=300)
plt.savefig('forge1.pdf',dpi=300)
plt.show()


def calculate_mse(original, reconstructed):
    """Calculate Mean Squared Error"""
    return np.mean((original - reconstructed) ** 2)

def calculate_rmse(original, reconstructed):
    """Calculate Root Mean Square Error"""
    return np.sqrt(calculate_mse(original, reconstructed))

def calculate_psnr(original, reconstructed):
    """Calculate Peak Signal-to-Noise Ratio"""
    # Add a small constant to avoid division by zero
    mse = calculate_mse(original, reconstructed)
    if mse == 0:
        return np.inf  # If MSE is zero, PSNR is infinite
    max_pixel = np.max(original)
    return 20 * np.log10(max_pixel / np.sqrt(mse))

def calculate_ssim(original, reconstructed):
    """Calculate Structural Similarity Index (SSIM)"""
    from skimage.metrics import structural_similarity as ssim
    ssim_value, _ = ssim(original, reconstructed, full=True)
    return ssim_value


# Generate synthetic data
original_data = data_dphase_bp
reconstructed_data = wav1
# Calculate metrics
mse_value = calculate_mse(original_data, reconstructed_data)
rmse_value = calculate_rmse(original_data, reconstructed_data)
psnr_value = calculate_psnr(original_data, reconstructed_data)
ssim_value = calculate_ssim(original_data, reconstructed_data)


#tensorflow
tf.placeholder() is not compatible with eager execution.
Solution: 
# import tensorflow as tf
import tensorflow.compat.v1 as tf 
tf.compat.v1.disable_eager_execution() 

#ML
import numpy as np;import keras;
inputs = np.random.random((32, 10, 8))
lstm = keras.layers.LSTM(4)
output = lstm(inputs)
print('output.shape',output.shape)

lstm = keras.layers.LSTM(4, return_sequences=True, return_state=True)
whole_seq_output, final_memory_state, final_carry_state = lstm(inputs)
print('whole_seq_output.shape',whole_seq_output.shape)
print('final_memory_state.shape',final_memory_state.shape)
print('final_carry_state.shape',final_carry_state.shape)

#github readme
<p align="center">
<img src='https://github.com/chenyk1990/gallery/blob/main/pywave/vel3d.png' alt='comp' width=960/>
<img src='https://github.com/chenyk1990/gallery/blob/main/pywave/data3d.png' alt='comp' width=640/>
</p>

#GIF/gif file from pngs
# more detailed example in pywave/demo/test_second.py
from PIL import Image

def create_gif(image_paths, output_gif_path, duration=500):
    """Creates a GIF from a list of PNG images."""

    images = [Image.open(image_path) for image_path in image_paths]

    images[0].save(
        output_gif_path,
        save_all=True,
        append_images=images[1:],
        optimize=False,
        duration=duration,
        loop=0  # 0 means infinite loop
    )

## create GIF/gif
inpath=fignames
outpath='wfd3ds.gif'
create_gif(inpath,outpath)
    

## play with the PIL, image to numpy coversion
	import PIL
	import numpy as np
	I = np.asarray(PIL.Image.open('dn.png'))
	im = PIL.Image.fromarray(np.uint8(I))
	im.save('dn2.png')


### EQCCT
def create_cct_modelP(inputs):

    inputs1 = convF1(inputs,   10, 11, 0.1)
    inputs1 = convF1(inputs1, 20, 11, 0.1)
    inputs1 = convF1(inputs1, 40, 11, 0.1)
    
    inputreshaped = layers.Reshape((6000,1,40))(inputs1)

    # Create patches.
    patches = Patches(patch_size)(inputreshaped)
    
    # Encode patches.
    encoded_patches = PatchEncoder(num_patches, projection_dim)(patches)
        
    # Calculate Stochastic Depth probabilities.
    dpr = [x for x in np.linspace(0, stochastic_depth_rate, transformer_layers)]

    # Create multiple layers of the Transformer block.
    for i in range(transformer_layers):
        # Layer normalization 1.
        x1 = layers.LayerNormalization(epsilon=1e-6)(encoded_patches)

        # Create a multi-head attention layer.
        attention_output = layers.MultiHeadAttention(
            num_heads=num_heads, key_dim=projection_dim, dropout=0.1
        )(x1, x1)

        # Skip connection 1.
        attention_output = StochasticDepth(dpr[i])(attention_output)
        x2 = layers.Add()([attention_output, encoded_patches])

        # Layer normalization 2.
        x3 = layers.LayerNormalization(epsilon=1e-6)(x2)

        # MLP.
        x3 = mlp(x3, hidden_units=transformer_units, dropout_rate=0.1)

        # Skip connection 2.
        x3 = StochasticDepth(dpr[i])(x3)
        encoded_patches = layers.Add()([x3, x2])

    # Apply sequence pooling.
    representation = layers.LayerNormalization(epsilon=1e-6)(encoded_patches)

    return representation


def create_cct_modelS(inputs):

    inputs1 = convF1(inputs,   10, 11, 0.1)
    inputs1 = convF1(inputs1, 20, 11, 0.1)
    inputs1 = convF1(inputs1, 40, 11, 0.1)
    
    inputreshaped = layers.Reshape((6000,1,40))(inputs1)

    # Create patches.
    patches = Patches(patch_size)(inputreshaped)
    
    # Encode patches.
    encoded_patches = PatchEncoder(num_patches, projection_dim)(patches)
        
    # Calculate Stochastic Depth probabilities.
    dpr = [x for x in np.linspace(0, stochastic_depth_rate, transformer_layers)]

    # Create multiple layers of the Transformer block.
    for i in range(transformer_layers):
        encoded_patches = convF1(encoded_patches, 40,11, 0.1)
        # Layer normalization 1.
        x1 = layers.LayerNormalization(epsilon=1e-6)(encoded_patches)

        # Create a multi-head attention layer.
        attention_output = layers.MultiHeadAttention(
            num_heads=num_heads, key_dim=projection_dim, dropout=0.1
        )(x1, x1)
        attention_output = convF1(attention_output, 40,11, 0.1)
    
        # Skip connection 1.
        attention_output = StochasticDepth(dpr[i])(attention_output)
        x2 = layers.Add()([attention_output, encoded_patches])

        # Layer normalization 2.
        x3 = layers.LayerNormalization(epsilon=1e-6)(x2)

        # MLP.
        x3 = mlp(x3, hidden_units=transformer_units, dropout_rate=0.1)

        # Skip connection 2.
        x3 = StochasticDepth(dpr[i])(x3)
        encoded_patches = layers.Add()([x3, x2])

    # Apply sequence pooling.
    representation = layers.LayerNormalization(epsilon=1e-6)(encoded_patches)

    return representation


#legend
h1=plt.plot(paths[0,:],paths[1,:],paths[2,:],'k',markersize=20)
h2=plt.plot(paths1[0,:],paths1[1,:],paths1[2,:],'g--',markersize=20)
# Plots endpoint
plt.legend([h1[0],h2[0]],['Ray 1', 'Ray 2'], loc='upper left')



