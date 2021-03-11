import pandas as pd  
import numpy as np  
import sys
import matplotlib.pyplot as plt  
import seaborn as sb 
from sklearn.model_selection import train_test_split 
from sklearn.linear_model import LinearRegression
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
print(df.describe())

x_cols = ["avg_exc", "avg_inh", "max_exc", "max_inh", "num_exc", "num_inh", "std_exc", "std_inh", "skew_exc", "skew_inh"]

for col in x_cols:
    plt.figure()
    sb.distplot(df[col])
    plt.show()
    
X = df[x_cols].values
y = df['FR'].values

vals = np.array(y).astype("int")

fractions = []
for i in range(0, int(np.max(y)) + 5):
    fractions.append(len(np.where(np.array(vals) == i)[0]) / len(y))

plt.figure()
plt.plot(fractions)
plt.show()

plt.figure()
plt.hist(np.log(np.array(y) + 1))
plt.show()

plt.figure()
sb.distplot(df['FR'])
plt.show()

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)

regressor = LinearRegression()
regressor.fit(X_train, y_train)

coefs = regressor.coef_
weighted_coefs = coefs.copy()
for i in range(len(coefs)):
    weighted_coefs[i] *= np.std(np.abs(df[x_cols[i]]))

#coeff_df = pd.DataFrame({"Coefs":regressor.coef_, "W_coefs": weighted_coefs}, x_cols)
coeff_df = pd.DataFrame(regressor.coef_, x_cols, columns=['Coefficients'])
print(coeff_df)

y_pred = regressor.predict(X_test)

plt.figure()
sb.distplot(df['FR'], label="trues")
sb.distplot(y_pred, label="preds")
plt.ylabel("Fraction of cells")
plt.xlabel("Firing rate")
plt.legend()
plt.show()

# comp = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
# print(comp.head(25))

# plt.figure()
# plt.hist(y_pred, label="preds", alpha = 0.5)
# plt.hist(y_test, label='trues', alpha = 0.5)
# plt.legend()
# plt.show()

# plt.figure()
# comp.head(25).plot(kind='bar',figsize=(10,8))
# plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
# plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
# plt.show()

print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))  
print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))  
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, y_pred)))
print("R2 Score:", metrics.r2_score(y_test, y_pred))
