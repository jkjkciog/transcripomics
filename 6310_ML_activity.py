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

Y=df.iloc['class']
print(Y.head())