import pandas as pd  
import numpy as np  
import sys
import matplotlib.pyplot as plt  
import seaborn as sb 
from sklearn.model_selection import train_test_split 
from sklearn.ensemble import RandomForestRegressor
from sklearn import metrics

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        inp = sys.argv[-1]
    else:
        raise Exception("no work" + str(sys.argv[-1]))

fname = str(inp)

#df = pd.read_csv("all_log_results.csv", index_col='gid')
#df = pd.read_csv("all_log_dend_results.csv", index_col='gid')
df = pd.read_csv(fname, index_col='gid')
df.shape

x_cols = ["avg_exc", "avg_inh", "max_exc", "max_inh", "num_exc", "num_inh", "std_exc", "std_inh", "skew_exc", "skew_inh"]

for col in x_cols:
    plt.figure()
    sb.distplot(df[col])
    plt.show()
    
X = df[x_cols].values
y = df['FR'].values

vals = np.array(y).astype("int")

# plt.figure()
# plt.hist(np.log(np.array(y) + 1))
# plt.show()

plt.figure()
sb.distplot(df['FR'])
plt.show()

plt.figure()
sb.distplot(df['FR'])
plt.xscale('log')
plt.show()

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)

# Instantiate model with 1000 decision trees
rf = RandomForestRegressor(n_estimators = 100, random_state = 42)
# Train the model on training data
rf.fit(X_train, y_train)

y_pred = rf.predict(X_test)

plt.figure()
sb.distplot(df['FR'], label="trues")
sb.distplot(y_pred, label="preds")
plt.ylabel("Fraction of cells")
plt.xlabel("Firing rate")
plt.legend()
plt.show()

print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))  
print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))  
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, y_pred)))
print("R2 Score:", metrics.r2_score(y_test, y_pred))

# Get numerical feature importances
importances = list(rf.feature_importances_)
# List of tuples with variable and importance
feature_importances = [(feature, round(importance, 2)) for feature, importance in zip(x_cols, importances)]
# Sort the feature importances by most important first
feature_importances = sorted(feature_importances, key = lambda x: x[1], reverse = True)
# Print out the feature and importances 
[print('Variable: {:20} Importance: {}'.format(*pair)) for pair in feature_importances];
