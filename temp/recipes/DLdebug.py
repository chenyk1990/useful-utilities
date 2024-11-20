#(python-version,package-version,date)

from keras.layers.core import Dense, Dropout, Activation, Flatten
-> from keras.layers import Dense, Dropout, Activation, Flatten
(p3.12,keras-3.6.0,202411)

from keras.layers.convolutional import Convolution2D, MaxPooling2D
-> from keras.layers import Convolution2D, MaxPooling2D
(p3.12,keras-3.6.0,202411)

from keras.utils import np_utils (np_utils.to_categorical) #this works in 
->from keras.utils import to_categorical
(p3.12,keras-3.6.0,202411)

model.add(Convolution2D(32, 1, 3, 3, border_mode='full'))
->model.add(Convolution2D(32, 3) or model.add(Convolution2D(32, 3, activation='relu'))
(p3.12,keras-3.6.0,202411)

model.add(MaxPooling2D(poolsize=(2, 2)))
->model.add(MaxPooling2D(pool_size=(2, 2)))
(p3.12,keras-3.6.0,202411)

model.add(Dense(1280,128))->model.add(Dense(128))
(p3.12,keras-3.6.0,202411)

model.fit(X_train, Y_train, batch_size=batch_size, nb_epoch=nb_epoch, show_accuracy=True, verbose=1, validation_data=(X_test, Y_test))
->model.compile(loss='categorical_crossentropy', optimizer='adadelta', metrics=['accuracy']) #show acc
model.fit(X_train, Y_train, batch_size=batch_size, epochs=nb_epoch, verbose=1, validation_data=(X_test, Y_test))
(p3.12,keras-3.6.0,202411)

score = model.evaluate(X_test, Y_test, show_accuracy=True, verbose=0)
->score = model.evaluate(X_test, Y_test, verbose=1)
(p3.12,keras-3.6.0,202411)

X_train = X_train.reshape(X_train.shape[0], 1, 28, 28) 		#torch NCHW
->X_train = X_train.reshape(X_train.shape[0], 28, 28, 1)	 #tf NHWC
(p3.12,keras-3.6.0,202411) 


# Equal commands
model.add(Convolution2D(32, 3));model.add(Activation('relu'))
=model.add(Convolution2D(32, 3,activation='relu'))
=model.add(Conv2D(32, 3,activation='relu'))