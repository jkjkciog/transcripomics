try:
    import numpy as np
    import pandas as pd
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import accuracy_score, classification_report
    from sklearn import svm

except: 
    print(f'Error importing libraries')

column_name= ["pregnancies", "glucose", "blood_pressure", "skin_thickness", "insulin", "bmi", "diabetes_pedigree_function", "age", "outcome"]

df = pd.read_csv("C:\\Users\\jaan\\Downloads\\6310_ML_activity.csv", names=column_name)


print(df.shape)
print(df.head())

X=df.iloc[:,0:8]
print(X.head())

Y = df['outcome']
print(Y.head())

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.25, random_state=4)
print("X_train.shape = ", X_train.shape)
print("Y_train.shape = ", Y_train.shape)
print("X_test.shape = ", X_test.shape)
print("Y_test.shape = ", Y_test.shape)

print(X_test.head())

#normailze features X
scaler = StandardScaler()
scaler.fit(X_train)
X_train = scaler.transform(X_train)
print(X_train[:5,:])
scaler.fit(X_test)
X_test = scaler.transform(X_test)

#train model
clf = svm.SVC(kernel='linear')
print(clf.fit(X_train, Y_train))

#Decision accuracy
y_prediction = clf.predict(X_test)
Y_train_prediction = clf.predict(X_train)
print("Accuracy Score = ", accuracy_score(Y_test, y_prediction))
print(y_prediction)
y_ac = accuracy_score(Y_train, Y_train_prediction)

#Compairing different models
for k in ('linear', 'poly', 'rbf', 'sigmoid'):
    clf = svm.SVC(kernel=k)
    clf.fit(X_train, Y_train)
    y_prediction = clf.predict(X_test)
    print("Accuracy Score for ", k, " = ", accuracy_score(Y_test, y_prediction))