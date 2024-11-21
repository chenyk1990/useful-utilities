#(python-version,package-version,date)
from keras.layers.core import Dense, Dropout, Activation, Flatten, Dense, RepeatVector
-> from keras.layers import Dense, Dropout, Activation, Flatten, Dense, RepeatVector
(p3.12,keras-3.6.0,tf-2.16.2,202411)

from keras.layers.convolutional import Convolution2D, MaxPooling2D
-> from keras.layers import Convolution2D, MaxPooling2D
(p3.12,keras-3.6.0,tf-2.16.2,202411)

from keras.layers.convolutional import Conv2D, Conv2DTranspose, Conv1D, ZeroPadding1D)#this works in tensorflow==2.8.0;keras==2.8.0,p3.7.16 (EQCCT env)
-> from keras.layers import Conv2D, Conv2DTranspose, Conv1D, ZeroPadding1D, Input (or keras.Input)
(p3.12,keras-3.6.0,tf-2.16.2,202411)

from keras.layers.recurrent import LSTM
-> from keras.layers import LSTM,RNN,GRU
(p3.12,keras-3.6.0,tf-2.16.2,202411)

from keras.layers import recurrent (recurrent.JZS1)
->from keras.layers import RNN or ->from keras.src.layers import RNN
(p3.12,keras-3.6.0,tf-2.16.2,202411)

from keras.models import slice_X
->from keras.engine.training import _slice_arrays
->from tensorflow.python.keras.engine import _slice_arrays
(not working)
#maybe:
def slice_X(din,start,end):
	return din[int(start):int(end),:]

import keras.engine.training
->import tensorflow.python.keras.engine
(not working)

from keras.utils import np_utils (np_utils.to_categorical)#this works in tensorflow==2.8.0;keras==2.8.0,p3.7.16 (EQCCT env)
->from keras.utils import to_categorical
(p3.12,keras-3.6.0,tf-2.16.2,202411)

from keras.datasets.data_utils import get_file
->from keras.utils import get_file
(p3.12,keras-3.6.0,tf-2.16.2,202411)

model.add(Convolution2D(32, 1, 3, 3, border_mode='full'))
->model.add(Convolution2D(32, 3) or model.add(Convolution2D(32, 3, activation='relu'))
(p3.12,keras-3.6.0,tf-2.16.2,202411)

model.add(Convolution2D(32, 3, activation='relu', input_shape=(28, 28)))
->model.add(Convolution2D(32, 3, activation='relu', input_shape=(28, 28, 1)))
->model.add(Input(shape=(28,28,1)));model.add(Convolution2D(32, 3, activation='relu'))
(p3.12,keras-3.6.0,tf-2.16.2,202411)

model.add(MaxPooling2D(poolsize=(2, 2)))
->model.add(MaxPooling2D(pool_size=(2, 2)))
(p3.12,keras-3.6.0,tf-2.16.2,202411)

model.add(Dense(1280,128))->model.add(Dense(128))
(p3.12,keras-3.6.0,tf-2.16.2,202411)

model.fit(X_train, Y_train, batch_size=batch_size, nb_epoch=nb_epoch, show_accuracy=True, verbose=1, validation_data=(X_test, Y_test))
->model.compile(loss='categorical_crossentropy', optimizer='adadelta', metrics=['accuracy']) #show acc
model.fit(X_train, Y_train, batch_size=batch_size, epochs=nb_epoch, verbose=1, validation_data=(X_test, Y_test))
(p3.12,keras-3.6.0,tf-2.16.2,202411)

model.add(LSTM(len(chars), 512, return_sequences=True))
->model.add(LSTM(512, return_sequences=True))
(p3.12,keras-3.6.0,tf-2.16.2,202411)

model.add(RNN(len(chars), HIDDEN_SIZE))
->model.add(RNN(HIDDEN_SIZE))
(p3.12,keras-3.6.0,tf-2.16.2,202411)

score = model.evaluate(X_test, Y_test, show_accuracy=True, verbose=0)
->score = model.evaluate(X_test, Y_test, verbose=1)
(p3.12,keras-3.6.0,tf-2.16.2,202411)

X_train = X_train.reshape(X_train.shape[0], 1, 28, 28) 		#torch NCHW
->X_train = X_train.reshape(X_train.shape[0], 28, 28, 1)	 #tf NHWC
(p3.12,keras-3.6.0,tf-2.16.2,202411) 

X = np.zeros((len(sentences), maxlen, len(chars)), dtype=np.bool)
-> X = np.zeros((len(sentences), maxlen, len(chars)), dtype=np.bool_)

a = np.log(a)/temperature;a = np.exp(a)/np.sum(np.exp(a));np.argmax(np.random.multinomial(1,a,1))
-> a = np.asarray(a).astype("float64");a = np.log(a)/temperature;a = np.exp(a)/np.sum(np.exp(a));np.argmax(np.random.multinomial(1,a,1))

xrange -> range
(p3.12)


model.predict_classes->model.predict

# Equal commands
model.add(Convolution2D(32, 3));model.add(Activation('relu'))
=model.add(Convolution2D(32, 3, activation='relu'))
=model.add(Convolution2D(32, 3, activation='relu', input_shape=(28, 28, 1)))
=model.add(Conv2D(32, 3,activation='relu'))



# install

pip install scikit-image  (skimage)
pip install tensorflow













