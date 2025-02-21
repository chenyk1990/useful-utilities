#(python-version,package-version,date)
from keras.layers.core import Dense, Dropout, Activation, Flatten, Dense, RepeatVector
-> from keras.layers import Dense, Dropout, Activation, Flatten, Dense, RepeatVector
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
from keras.layers.convolutional import Convolution2D, MaxPooling2D
-> from keras.layers import Convolution2D, MaxPooling2D
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
from keras.layers.convolutional import Conv2D, Conv2DTranspose, Conv1D, ZeroPadding1D)#this works in tensorflow==2.8.0;keras==2.8.0,p3.7.16 (EQCCT env)
-> from keras.layers import Conv2D, Conv2DTranspose, Conv1D, ZeroPadding1D, Input (or keras.Input)
concatenate, LeakyReLU, BatchNormalization
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
from keras.layers.recurrent import LSTM
-> from keras.layers import LSTM,RNN,GRU
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
from keras.layers.advanced_activations import PReLU
->from keras.layers import PReLU
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
from keras.layers import recurrent (recurrent.JZS1)
->from keras.layers import RNN or ->from keras.src.layers import RNN
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
from keras.models import slice_X
->from keras.engine.training import _slice_arrays
->from tensorflow.python.keras.engine import _slice_arrays
(not working)
#maybe:
def slice_X(din,start,end):
	return din[int(start):int(end),:]
###############################################################################################
import keras.engine.training
->import tensorflow.python.keras.engine
(not working)
###############################################################################################
from keras.utils import np_utils (np_utils.to_categorical)#this works in tensorflow==2.8.0;keras==2.8.0,p3.7.16 (EQCCT env)
->from keras.utils import to_categorical
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
from keras.datasets.data_utils import get_file
->from keras.utils import get_file
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
model.add(PReLU((512,)))->model.add(PReLU())
model.add(BatchNormalization((512,))) -> model.add(BatchNormalization())
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
model.add(Convolution2D(32, 1, 3, 3, border_mode='full'))
->model.add(Convolution2D(32, 3) or model.add(Convolution2D(32, 3, activation='relu'))
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
model.add(Convolution2D(32, 3, activation='relu', input_shape=(28, 28)))
->model.add(Convolution2D(32, 3, activation='relu', input_shape=(28, 28, 1)))
->model.add(Input(shape=(28,28,1)));model.add(Convolution2D(32, 3, activation='relu'))
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
model.add(MaxPooling2D(poolsize=(2, 2)))
->model.add(MaxPooling2D(pool_size=(2, 2)))
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
model.add(Dense(1280,128))->model.add(Dense(128))
(p3.12,keras-3.6.0,tf-2.16.2,202411)

model.add(Dense(dims, 512, init='glorot_uniform'))
->model.add(Dense(512))
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
model.fit(X_train, Y_train, batch_size=batch_size, nb_epoch=nb_epoch, show_accuracy=True, verbose=1, validation_data=(X_test, Y_test))
->model.compile(loss='categorical_crossentropy', optimizer='adadelta', metrics=['accuracy']) #show acc
model.fit(X_train, Y_train, batch_size=batch_size, epochs=nb_epoch, verbose=1, validation_data=(X_test, Y_test))
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
def R2_score(v_true, v_pred):
    from keras import backend as K
    ssres = K.sum(K.square(v_true - v_pred))
    sstot = K.sum(K.square(v_true - K.mean(v_true)))
    return 1 - ssres / sstot #or from sklearn.metrics import r2_score (coefficient of determination) regression score function.)
c = keras.optimizers.adam(lr = lr_val);cnn_model.compile(optimizer=c, loss='mse', metrics=[R2_score])
->cnn_model.compile(loss='mse', optimizer='adam', metrics=[R2_score])
(EQCCT env)
###############################################################################################
model.add(LSTM(len(chars), 512, return_sequences=True))
->model.add(LSTM(512, return_sequences=True))
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
model.add(RNN(len(chars), HIDDEN_SIZE))
->model.add(RNN(HIDDEN_SIZE))
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
score = model.evaluate(X_test, Y_test, show_accuracy=True, verbose=0)
->score = model.evaluate(X_test, Y_test, verbose=1)
(p3.12,keras-3.6.0,tf-2.16.2,202411)
###############################################################################################
X_train = X_train.reshape(X_train.shape[0], 1, 28, 28) 		#torch NCHW
->X_train = X_train.reshape(X_train.shape[0], 28, 28, 1)	 #tf NHWC
(p3.12,keras-3.6.0,tf-2.16.2,202411) 
###############################################################################################
X = np.zeros((len(sentences), maxlen, len(chars)), dtype=np.bool)
-> X = np.zeros((len(sentences), maxlen, len(chars)), dtype=np.bool_)
###############################################################################################
a = np.log(a)/temperature;a = np.exp(a)/np.sum(np.exp(a));np.argmax(np.random.multinomial(1,a,1))
-> a = np.asarray(a).astype("float64");a = np.log(a)/temperature;a = np.exp(a)/np.sum(np.exp(a));np.argmax(np.random.multinomial(1,a,1))
###############################################################################################
xrange -> range
(p3.12)
###############################################################################################
#downgrade
from keras.src import ops
ops.matmul -> import tensorflow as tf; tf.matmul
(keras-3.6.0->keras-2.11.0) 
###############################################################################################
model.predict_classes->model.predict
model.predict_proba -> model.predict
###############################################################################################
# Equal commands
model.add(Convolution2D(32, 3));model.add(Activation('relu'))
=model.add(Convolution2D(32, 3, activation='relu'))
=model.add(Convolution2D(32, 3, activation='relu', input_shape=(28, 28, 1)))
=model.add(Conv2D(32, 3,activation='relu'))




#torch
from torch.utils.tensorboard import SummaryWriter -> 
(Using: conda install -c conda-forge tensorboard)














# install
pip install scikit-image  (skimage)
pip install tensorflow
pip install scikit-learn
pip install torchvision
pip install h5py

#common
from sklearn.preprocessing import LabelEncoder,StandardScaler
from keras.models import Sequential
from keras.layers import Activation, Dense, RepeatVector, Input
from keras.layers import RNN
from keras.src import ops
from sklearn.utils import shuffle
import numpy as np
model = Sequential()
model.compile(loss='categorical_crossentropy', optimizer='adadelta', metrics=['accuracy'])
help(model.compile)
###############################################################################################
#generate simple data
import numpy as np
from pyseistr import gensyn
data,noisy=gensyn(noise=True);[n1,n2]=data.shape;
# import matplotlib.pyplot as plt;
# plt.subplot(1,2,1);plt.imshow(data,clim=[-0.2,0.2],aspect='auto');plt.xlabel('Trace');plt.ylabel('Time sample');
# plt.subplot(1,2,2);plt.imshow(noisy,clim=[-0.2,0.2],aspect='auto');plt.xlabel('Trace');
# plt.show();
from pyseistr import patch2d,patch2d_inv,snr
X=patch2d(data,l1=16,l2=16,s1=8,s2=8);
Xnoisy=patch2d(noisy,l1=16,l2=16,s1=8,s2=8);
from keras import layers
from keras.models import Model
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau
input = layers.Input(shape=(256,))  #or from tensorflow.keras.layers import Input, Dense
x = layers.Dense(64, activation="relu", name="layer1")(input)  
x = layers.Dense(256, activation="linear", name="layer2")(x)  
# Autoencoder
model = Model(input,x)
model.summary()
# Compile
model.compile(optimizer="adam", loss='mse')
import datetime
today=datetime.date.today()
weightname='best_model_%s.weights.h5'%str(today)
checkpoint = ModelCheckpoint(filepath=weightname,monitor='val_loss',mode = 'min',verbose=1,save_weights_only=True,save_best_only=True)
lr_reducer = ReduceLROnPlateau(factor=0.1,cooldown=0,patience=50,min_lr=0.5e-6,monitor='val_loss',mode = 'min',verbose= 1)
model.fit(Xnoisy,Xnoisy,batch_size=128,verbose=1,epochs=20,callbacks=[checkpoint,lr_reducer],validation_split=0.2)
###############################################################################################
import os
os.environ["KERAS_BACKEND"] = "tensorflow"


###############################################################################################
#torch #tensorflow #Height#Width
Pytorch input dimension:
For a conv2D, input should be in (N, C, H, W) format. N is the number of samples/batch_size. C is the channels. H and W are height and width resp.
For conv1D, input should be (N,C,L) see documentation at
TF ipnut dimension
conv2D: (N,H,W,C) (a little different)

#torch#tensorflow#Height
Input and Output (Weight and height) in convolutional layer
W2=(W1+2P-K)/S+1
H2=(H1+2P-K)/S+1
#number of parameters
CONV layer: This is where CNN learns, so certainly we’ll have weight matrices. To calculate the learnable parameters here, all we have to do is just multiply the by the shape of width m, height n, previous layer’s filters d and account for all such filters k in the current layer. Don’t forget the bias term for each of the filter. Number of parameters in a CONV layer would be : ((m * n * d)+1)* k), added 1 because of the bias term for each filter. The same expression can be written as follows: ((shape of width of the filter * shape of height of the filter * number of filters in the previous layer+1)*number of filters). Where the term “filter” refer to the number of filters in the current layer.
###############################################################################################





from keras.optimizers import Adam

