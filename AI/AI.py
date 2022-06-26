# sourcery skip: avoid-builtin-shadow
import numpy as np
import pandas
from imblearn.over_sampling import SMOTE
from sklearn.tree import ExtraTreeClassifier

import tensorflow as tf
import datetime

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from tensorflow.python.keras.layers import Dense, Dropout
from tensorflow.python.keras.models import Sequential

from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.naive_bayes import MultinomialNB
from sklearn.neighbors import KNeighborsClassifier
from tensorboard.plugins.hparams import api as hp


data = []
with open("prepared_data\\train_data.csv") as f:
    lis = [line.split(",") for line in f]
    for i in lis:
        data.append([float(idx) for idx in i])

dataframe = pandas.DataFrame(data)
dataset = dataframe.values
eff = dataset[:, 0:1500].astype(float)

data = []
with open("prepared_data\\train_data_neg.csv") as f:
    for str in f:
        str = str.replace('"', '')
        str = str[:-2]
        data.append([float(idx) for idx in str.split(',')])

dataframe = pandas.DataFrame(data)
dataset = dataframe.values
noneff = dataset[:, 0:1500].astype(float)


data = []
with open("prepared_data\\train_data_all.csv") as f:
    for str in f:
        str = str.replace('"', '')
        str = str[:-2]
        data.append([float(idx) for idx in str.split(',')])


dataframe = pandas.DataFrame(data)
dataset = dataframe.values
featurevector = dataset[:, 0:1500].astype(float)

print("1")

a1 = eff.shape
a2 = noneff.shape
Y = np.ones((a1[0] + a2[0], 1))

for i in range(a1[0]):
    Y[i, 0] = 0

print('Resampling the unbalanced data...')
X_resampled, Y_resampled = SMOTE().fit_resample(featurevector, Y)

scaler = MinMaxScaler().fit(X_resampled)
X_resampled = scaler.transform(X_resampled)

featurevector = scaler.transform(featurevector)
newshape = X_resampled.shape

X_train, X_test, y_train, y_test = train_test_split(X_resampled, Y_resampled, test_size=0.15)
y_t = y_train
y_te = y_test
y_train = np.ones((len(y_t), 2))
y_test = np.ones((len(y_te), 2))

for i in range(len(y_t)):
    if y_t[i] == 0:
        y_train[i][1] = 0
    if y_t[i] == 1:
        y_train[i][0] = 0

for i in range(len(y_te)):
    if y_te[i] == 0:
        y_test[i][1] = 0
    if y_te[i] == 1:
        y_test[i][0] = 0        

print(2)

print("Training Artificial Neural Network...")
model = Sequential()
model.add(Dense(newshape[1] + 1, activation='relu', input_shape=(newshape[1],)))
model.add(Dense(750, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(350, activation='relu'))
model.add(Dense(125, activation='relu'))
model.add(Dense(75, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(30, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(2, activation='sigmoid'))
model.compile(loss='binary_crossentropy',optimizer='adam',metrics=['binary_accuracy'])


log_dir="logs/fit/" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=log_dir, histogram_freq=1)
hparams_callback = hp.KerasCallback(log_dir, {'num_relu_units': 512,'dropout': 0.2})

model.fit(X_train, y_train, epochs=20, batch_size=25, verbose=1,callbacks=[tensorboard_callback, hparams_callback])
score = model.evaluate(X_test, y_test, verbose=1)
ANN = model.predict(X_test)
ANN = model.predict(featurevector)
model.summary()



y_train = []
y_test = []
y_train = y_t
y_test = y_te

print("Training Support Vector Machine...")
clf1 = svm.SVC(decision_function_shape='ovr', kernel='linear', max_iter=1000, probability=True)
clf1.fit(X_train, y_train)
y_pred = clf1.predict(X_test)
results = cross_val_score(clf1, X_test, y_test, cv=10)
SVM = clf1.predict_proba(featurevector)

print("Training k-Nearest Neighbor ...")
neigh = KNeighborsClassifier(n_neighbors=10)
neigh.fit(X_train, y_train)
results = cross_val_score(neigh, X_test, y_test, cv=10)
y_pred = neigh.predict(X_test)
KNN = neigh.predict_proba(featurevector)


print("Training Naive Bayes...")
clf = MultinomialNB()
clf.fit(X_train, y_train)
results = cross_val_score(clf, X_test, y_test, cv=10)
y_pred = clf.predict(X_test)
DT = clf.predict_proba(featurevector)

print("Training Random Forest...")
rf = RandomForestClassifier(random_state=0, min_samples_leaf=100)
rf.fit(X_train, y_train)
results = cross_val_score(rf, X_test, y_test, cv=10)
y_pred = rf.predict(X_test)
RF = clf.predict_proba(featurevector)

vote_result = [[] for _ in range(len(SVM))]
for i in range(len(ANN)):
    vote_result[i].append(ANN[i][0])
    vote_result[i].append(SVM[i][0])
    vote_result[i].append(KNN[i][0])
    vote_result[i].append(DT[i][0])
    vote_result[i].append(RF[i][0])

for i in range(len(ANN)):
    res = np.mean(vote_result[i])
    if res >= 0.8:
        print('Sequence ', i + 1, ' is a probable Type 6 Effector')
    else:
        print('Sequence ', i + 1, ' is not a Type 6 Effector')
