import numpy as np
import pandas as pd
# import analyze_spiking
import scipy.signal as s
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.model_selection import train_test_split
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import my_plotting

#Bins spiking by ms and filters for a smooth gamma oscillation.
def generate_spike_gamma(file, time):
    data = my_plotting.load_dataset(file)
    timestamps = np.array(data['timestamps'])

    bins = [len(np.where(np.trunc(timestamps) == i)[0]) for i in range(time)]

    b, a = s.butter(2, [50, 80], btype = "bandpass",fs = 1000)

    gamma = s.filtfilt(b, a, bins)

    return gamma

syn_info = pd.read_csv("exc_syn_info.csv", index_col="node_id")
spike_file = "gamma_results/spikes.h5"
time = 120000
inh_file = 'gamma_results/inh_stim_spikes.h5'
exc_file = 'exc_stim_spikes.h5'

def make_data():
    exc_inp = my_plotting.load_dataset(exc_file)
    exc_ts = np.array(exc_inp['timestamps'])
    exc_ids = np.array(exc_inp['node_ids'])

    spikes = np.array(my_plotting.load_dataset(spike_file)['timestamps'])
    gamma = generate_spike_gamma(inh_file, time)
    troughs = s.find_peaks(-gamma)[0]

    #import pdb; pdb.set_trace()
    exc_phases = np.zeros(len(exc_ts))
    for i in range(len(troughs) - 1):
        start = troughs[i]
        stop = troughs[i + 1]
        length = stop - start

        aps = np.where((exc_ts >= start) & (exc_ts < stop))[0]
        times = exc_ts[aps]
        times = times - start
        times = times / length
        times = np.abs(times - 0.5)
        exc_phases[aps] = times

    weights = np.array([syn_info.loc[id]["weight"] for id in exc_ids])
    distance = np.array([syn_info.loc[id]["distance"] for id in exc_ids])
    is_basal = np.array([syn_info.loc[id]["is_basal"] for id in exc_ids])

    ids_s = []
    last = -1
    for t in spikes:
        if last == -1 or t > last + 100:
            #ids = np.where((exc_ts >= t - 20) & (exc_ts < t))[0]
            #nums.append(len(ids))
            #spike_follows[ids] = True
            #print(t)
            ids_s.append(t - 20)
        last = t
    
    ids = []
    import random
    random.seed(42)
    while len(ids) < len(ids_s) :#* 10:
        id = np.random.randint(0, time)
        if id not in ids:
            if len(np.where((spikes >= id) & (spikes < id + 40))[0]) == 0 and len(np.where((exc_ts >= id) & (exc_ts < id + 20))[0]) != 0:
                ids.append(id)
        if (len(ids) >= len(ids_s)):
            break
    
    ns = []#num spikes
    ws = []#mean weight
    ds = []#mean distance
    ibs = []#prop basal
    ps = []#mean phase
    sf = []#before spike
    for id in ids_s:
        inps = np.where((exc_ts >= id) & (exc_ts < id + 20))[0]
        ns.append(len(inps))
        ws.append(np.mean(weights[inps]))
        ds.append(np.mean(distance[inps]))

        is_basals = is_basal[inps]
        ibs.append(len(np.where(is_basals == True)[0]) / len(is_basals))
        ps.append(np.mean(exc_phases[inps]))
        #phases = exc_phases[inps]

        sf.append(True)

    for id in ids:
        inps = np.where((exc_ts >= id) & (exc_ts < id + 20))[0]
        ns.append(len(inps))
        ws.append(np.mean(weights[inps]))
        ds.append(np.mean(distance[inps]))

        is_basals = is_basal[inps]
        ibs.append(len(np.where(is_basals == True)[0]) / len(is_basals))
        ps.append(np.mean(exc_phases[inps]))
        sf.append(False)


    #import pdb; pdb.set_trace()

    df = pd.DataFrame()
    # df['weight'] = weights
    # df['distance'] = distance
    # df['is_basal'] = is_basal
    # df['phase'] = exc_phases
    # df['spike'] = spike_follows
    # df.to_csv("logistic_data.csv", index=False)
    df['num_spikes'] = ns
    df['weight'] = ws
    df['distance'] = ds
    df['is_basal'] = ibs
    df['phase'] = ps
    df['spike'] = sf
    df.to_csv("logistic_data2.csv", index=False)

# make_data()
# import pdb; pdb.set_trace()

df = pd.read_csv("logistic_data2.csv")
print(df.describe())

#X = df[["weight", "distance", "is_basal", "phase", "num_spikes"]]
X = df[["distance", "num_spikes", "phase"]]
y = df['spike']

# from imblearn.over_sampling import SMOTE
# os = SMOTE(random_state=0)
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=0)
# columns = X_train.columns
# X_train,y_train=os.fit_sample(X_train, y_train)

import statsmodels.api as sm
logit_model=sm.Logit(y,X.astype(float))
result=logit_model.fit()
print(result.summary2())

#import pdb; pdb.set_trace()
# we can Check the numbers of our data
# print("length of oversampled data is ",len(os_data_X))
# print("Number True:",len(os_data_y[os_data_y['y']==True]))
# print("Number False:",len(os_data_y[os_data_y['y']==False]))
# print("Proportion of no subscription data in oversampled data is ",len(os_data_y[os_data_y['y']==False])/len(os_data_X))
# print("Proportion of subscription data in oversampled data is ",len(os_data_y[os_data_y['y']==True])/len(os_data_X))
#import pdb; pdb.set_trace()
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=0)
logreg = LogisticRegression()
logreg.fit(X_train, y_train)

y_pred = logreg.predict(X_test)
print('Accuracy of logistic regression classifier on test set: {:.2f}'.format(logreg.score(X_test, y_test)))

from sklearn.metrics import confusion_matrix
confusion_matrix = confusion_matrix(y_test, y_pred)
print(confusion_matrix)

from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
import matplotlib.pyplot as plt
logit_roc_auc = roc_auc_score(y_test, logreg.predict(X_test))
fpr, tpr, thresholds = roc_curve(y_test, logreg.predict_proba(X_test)[:,1])
plt.figure()
plt.plot(fpr, tpr, label='Logistic Regression (area = %0.2f)' % logit_roc_auc)
plt.plot([0, 1], [0, 1],'r--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.savefig('Log_ROC')
plt.show()

import pdb; pdb.set_trace()