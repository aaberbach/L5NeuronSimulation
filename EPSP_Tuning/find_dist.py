#This code is used to create a distribution of uEPSPs that matches Song et al. 2005

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb

def lognormal(m, s):
        mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
        std = np.sqrt(np.log((s/m)**2 + 1))
        return max(np.random.lognormal(mean, std, 1), 0.0000000001)

m = 0.401
#s = 0.31
s=0.31

#import pdb; pdb.set_trace()
arr = [float(lognormal(m, s)) for i in range(50000)]
df = pd.DataFrame()
df['dist'] = arr

#import pdb; pdb.set_trace()
sb.distplot(df['dist'], hist=True)
plt.show()

sb.distplot(df['dist'], hist=True)
plt.xlabel('uEPSP (mv)')
plt.xscale('log')
plt.show()