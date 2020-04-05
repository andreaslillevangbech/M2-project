import numpy as np
import pandas as pd
from keras import models
from keras import layers, optimizers
from keras.callbacks import EarlyStopping, ModelCheckpoint
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt


df = pd.read_csv("pydf.csv")

print(df.head(10))

train_T = round(max(df['time']) / 2)
val_T = round(max(df['time']))*3/4
train = df[df['time'] <= train_T]
val = df[(train_T < df['time']) & (df['time'] <= val_T)]
test = df[df['time'] > val_T]
train.drop(['time', 'id'], axis = 1, inplace = True)
val.drop(['time', 'id'], axis = 1, inplace = True)
test.drop(['time', 'id'], axis = 1, inplace = True)

# Split into y and x
x_train = train.loc[:,'X1':]
y_train = train.loc[:,'Y']
x_val = val.loc[:,'X1':]
y_val = val.loc[:,'Y']
x_test = test.loc[:,'X1':]
y_test = test.loc[:,'Y']

# Standardize the inputs
scaler = StandardScaler().fit(x_train)
x_train = scaler.transform(x_train)
x_test = scaler.transform(x_test)


np.random.seed(0)

number_of_features = len(x_train[1]) 

# Start neural network
network = models.Sequential()

# Add fully connected layer with a ReLU activation function
network.add(layers.Dense(units=16, activation='relu', input_shape=(number_of_features,)))
network.add(layers.normalization.BatchNormalization())

# Add fully connected layer with a ReLU activation function
network.add(layers.Dense(units=8, activation='relu'))
network.add(layers.normalization.BatchNormalization())

# Add fully connected layer with a sigmoid activation function
network.add(layers.Dense(units=1, activation='linear'))

# Compile neural network
sgd = optimizers.SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
network.compile(loss='mse', # RSS
                optimizer='rmsprop', 
                metrics=['mae']) # Stochastic gradient descent

# Set callback functions to early stop training and save the best model so far
callbacks = [EarlyStopping(monitor='val_loss', patience=2),
             ModelCheckpoint(filepath='best_model.h5', monitor='val_loss', save_best_only=True)]

# Train neural network
history = network.fit(x_train, # Features
                      y_train, # Target vector
                      epochs=100, # Number of epochs
                      callbacks=callbacks, # Early stopping
                      verbose=0, # Print description after each epoch
                      batch_size=32, # Number of observations per batch
                      validation_data=(x_val, y_val)) # Data for evaluation

# Plotting the training progress
# Plot training & validation loss values
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('Model loss')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.legend(['Train', 'Test'], loc='upper left')
plt.show()


# Predict on test set
prediction = network.predict(x_test, batch_size=128)
mse = np.mean(np.square(prediction[:,0] - y_test))
total = np.mean(np.square(y_test - np.mean(y_train)))
R2 = (1 - mse/total)*100
print(R2)



np.savetxt('pred.NN3', prediction )

#mse.NN2 = np.mean((pred.NN2 - y_test)^2)

#print(test.mse)
