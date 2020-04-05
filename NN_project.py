import numpy as np
import pandas as pd
from keras import models
from keras import layers
from keras.callbacks import EarlyStopping, ModelCheckpoint

df = pd.read_csv("for_python.csv")

print(df.head(10))

train_T = round(max(df['time']) / 2)
val_T = round(max(df['time']))*3/4
train = df[df['time'] <= train_T]
val = df[(train_T < df['time']) & (df['time'] <= val_T)]
test = df[df['time'] > val_T]
train.drop(['time', 'id'], axis = 1, inplace = True)
val.drop(['time', 'id'], axis = 1, inplace = True)
test.drop(['time', 'id'], axis = 1, inplace = True)

np.random.seed(0)

number_of_features = len(train.iloc[1]) - 1 

# Start neural network
network = models.Sequential()

# Add fully connected layer with a ReLU activation function
network.add(layers.Dense(units=32, activation='relu', input_shape=(number_of_features,)))

# Add fully connected layer with a ReLU activation function
network.add(layers.Dense(units=16, activation='relu'))

# Add fully connected layer with a ReLU activation function
network.add(layers.Dense(units=8, activation='relu'))

# Add fully connected layer with a sigmoid activation function
network.add(layers.Dense(units=1, activation='linear')) 

# Compile neural network
sgd = optimizers.SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
network.compile(loss='mean_squared_error', # Cross-entropy
                optimizer='sgd') # Stochastic gradient descent

# Set callback functions to early stop training and save the best model so far
callbacks = [EarlyStopping(monitor='val_loss', patience=2),
             ModelCheckpoint(filepath='best_model.h5', monitor='val_loss', save_best_only=True)]

# Train neural network
history = network.fit(train_features, # Features
                      train_target, # Target vector
                      epochs=20, # Number of epochs
                      callbacks=callbacks, # Early stopping
                      verbose=0, # Print description after each epoch
                      batch_size=100, # Number of observations per batch
                      validation_data=(test_features, test_target)) # Data for evaluation