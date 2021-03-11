import numpy as np
import warnings
import scipy.signal as ss
import scipy
import scipy.stats as st
import matplotlib.pyplot as plt
import h5py
from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
import pandas as pd
from scipy.fft import fft
import matplotlib
import statsmodels.api as sm


# Create models from data
def best_fit_distribution(data, bins=200, ax=None):
    """Model data by finding best fit distribution to data"""
    # Get histogram of original data
    y, x = np.histogram(data, bins=bins, density=True)
    x = (x + np.roll(x, -1))[:-1] / 2.0

    # Distributions to check
    DISTRIBUTIONS = [        
        st.alpha,st.anglit,st.arcsine,st.beta,st.betaprime,st.bradford,st.burr,st.cauchy,st.chi,st.chi2,st.cosine,
        st.dgamma,st.dweibull,st.erlang,st.expon,st.exponnorm,st.exponweib,st.exponpow,st.f,st.fatiguelife,st.fisk,
        st.foldcauchy,st.foldnorm,st.frechet_r,st.frechet_l,st.genlogistic,st.genpareto,st.gennorm,st.genexpon,
        st.genextreme,st.gausshyper,st.gamma,st.gengamma,st.genhalflogistic,st.gilbrat,st.gompertz,st.gumbel_r,
        st.gumbel_l,st.halfcauchy,st.halflogistic,st.halfnorm,st.halfgennorm,st.hypsecant,st.invgamma,st.invgauss,
        st.invweibull,st.johnsonsb,st.johnsonsu,st.ksone,st.kstwobign,st.laplace,st.levy,st.levy_l,st.levy_stable,
        st.logistic,st.loggamma,st.loglaplace,st.lognorm,st.lomax,st.maxwell,st.mielke,st.nakagami,st.ncx2,st.ncf,
        st.nct,st.norm,st.pareto,st.pearson3,st.powerlaw,st.powerlognorm,st.powernorm,st.rdist,st.reciprocal,
        st.rayleigh,st.rice,st.recipinvgauss,st.semicircular,st.t,st.triang,st.truncexpon,st.truncnorm,st.tukeylambda,
        st.uniform,st.vonmises,st.vonmises_line,st.wald,st.weibull_min,st.weibull_max,st.wrapcauchy
    ]

    # Best holders
    best_distribution = st.norm
    best_params = (0.0, 1.0)
    best_sse = np.inf

    # Estimate distribution parameters from data
    for distribution in DISTRIBUTIONS:

        # Try to fit the distribution
        try:
            # Ignore warnings from data that can't be fit
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')

                # fit dist to data
                params = distribution.fit(data)

                # Separate parts of parameters
                arg = params[:-2]
                loc = params[-2]
                scale = params[-1]

                # Calculate fitted PDF and error with fit in distribution
                pdf = distribution.pdf(x, loc=loc, scale=scale, *arg)
                sse = np.sum(np.power(y - pdf, 2.0))

                # if axis pass in add to plot
                try:
                    if ax:
                        pd.Series(pdf, x).plot(ax=ax)
                    end
                except Exception:
                    pass

                # identify if this distribution is better
                if best_sse > sse > 0:
                    best_distribution = distribution
                    best_params = params
                    best_sse = sse

        except Exception:
            pass

    return (best_distribution.name, best_params)

def make_cdf(dist, params, size=10000):
    """Generate distributions's Probability Distribution Function """

    # Separate parts of parameters
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Get sane start and end points of distribution
    start = dist.ppf(0.01, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.01, loc=loc, scale=scale)
    end = dist.ppf(0.99, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.99, loc=loc, scale=scale)

    # Build PDF and turn into pandas Series
    x = np.linspace(start, end, size)
    y = dist.cdf(x, loc=loc, scale=scale, *arg)
    cdf = pd.Series(y, x)

    return cdf

def make_pdf(dist, params, size=10000):
    """Generate distributions's Probability Distribution Function """

    # Separate parts of parameters
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Get sane start and end points of distribution
    start = dist.ppf(0.02, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.01, loc=loc, scale=scale)
    end = dist.ppf(0.999, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.99, loc=loc, scale=scale)
    # Build PDF and turn into pandas Series
    x = np.linspace(start, end, size)
    y = dist.pdf(x, loc=loc, scale=scale, *arg)
    pdf = pd.Series(y, x)

    return pdf

def zscore(x):
    return (x-np.mean(x))/np.std(x)

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def make_noise(num_traces=100,num_samples=4999, mean_fr = 1, std_fr = 1):
    num_samples = num_samples+2000
    #fv = np.linspace(0, 1, 40);                                # Normalised Frequencies
    #a = 1/(1+2*fv);                                      # Amplitudes Of '1/f'
    #b = ss.firls(43, fv, a);                                   # Filter Numerator Coefficients
    B = [0.049922035, -0.095993537, 0.050612699, -0.004408786];
    A = [1, -2.494956002,   2.017265875,  -0.522189400];
    invfn = np.zeros((num_traces,num_samples))
    for i in np.arange(0,num_traces):
        wn = np.random.normal(loc=1,
 			   scale=0.5,size=num_samples)
        invfn[i,:] = zscore(ss.lfilter(B, A, wn));                             # Create '1/f' Noise
    return invfn[:,2000:]

def make_spikes(numUnits=100,rateProf=None):
    rateProf=rateProf
    rateProf[rateProf<0] = 0
    #meanRates = np.mean(rateProf)
    normRateProf= rateProf

    df_full = pd.DataFrame({'timestamps':[],'node_ids':[]})

    for i in np.arange(0,numUnits):
        rate_temp=[];simSpks_temp=[];
        rate_temp = normRateProf + np.exp(st.levy_stable.rvs(alpha=1.37, beta=-1.00, loc=0.92, scale=0.44, size=1))
        numbPoints = scipy.stats.poisson(rate_temp/1000).rvs()#Poisson number of points

        simSpks=np.where(numbPoints>0)[0]

        timestamps=simSpks
        node_ids = np.tile(i,simSpks.shape[0])
        df = pd.DataFrame({'timestamps':timestamps, 'node_ids':node_ids})
        df_full = df_full.append(df)

    return df_full

#### Load Experimental Data ####
data = pd.read_csv('L56FR_Drew.csv',header=None)
data = [np.log(i[0]) for i in data.values]
data = pd.Series(np.array(data).astype(float))


# Plot for comparison
#plt.figure(figsize=(12,8))
#ax = data.plot(kind='hist', bins=40, normed=True, alpha=0.5)
# Save plot limits
#dataYLim = ax.get_ylim()

# Find best fit distribution
#best_fit_name, best_fit_params = best_fit_distribution(data, 40, ax)
#best_dist = getattr(st, best_fit_name)

#print(best_fit_name, best_fit_params)
# Update plots
#ax.set_ylim(dataYLim)
#ax.set_title('{}{}'.format(best_fit_name,best_fit_params))
#ax.set_xlabel(u'Firing Rate (Hz)')
#ax.set_ylabel('Frequency')


# Make PDF with best params 
pdf = make_pdf(st.levy_stable, [1.37,-1.0,0.92,0.44])
cdf = make_cdf(st.levy_stable, [1.37,-1.0,0.92,0.44])

sim_data = st.levy_stable.rvs(alpha=1.37,beta=-1.0,loc=0.92,scale=0.44,size=1000)

# Display
#plt.figure()
#h = np.histogram(data,bins=40)
#plt.plot(h[1][1:],h[0]/np.sum(h[0]))
#plt.plot(pdf, lw=2, label='PDF',color='r')


#plt.figure()
#h = np.histogram(data,bins=40, density=True)
#plt.plot(h[1][1:],np.cumsum(h[0]/np.sum(h[0])))
#plt.plot(cdf, lw=2, label='CDF',color='r')

#param_names = (best_dist.shapes + ', loc, scale').split(', ') if best_dist.shapes else ['loc', 'scale']
#param_str = ', '.join(['{}={:0.2f}'.format(k,v) for k,v in zip(param_names, best_fit_params)])
#dist_str = '{}({})'.format(best_fit_name, param_str)
#ax.set_title('{}{}{}'.format(param_names, param_str,best_fit_name))

################################


z = make_noise(num_samples=20000-1,num_traces=100)


plt.figure()
plt.subplot(3,1,1)
plt.plot(z[0,:],linewidth=1)
plt.title('firing rate trace 1')

plt.subplot(3,1,2)
plt.hist(z[0,:],bins=100)

plt.subplot(3,1,3)
dt = 1/1000
from matplotlib.mlab import psd
psds = np.zeros((1,129))
for i in np.arange(0,1):
    psds[i,:], psd_frequencies = psd((z[i,:]-np.mean(z[i,:]))/np.std(z[i,:]), NFFT=256, Fs=1/dt)
plt.loglog(psd_frequencies, np.mean(psds,0))


################3

#df_full = make_spikes(numUnits=10000,rateProfMat=z)
#df_full2 = make_spikes(numUnits=1000,rateProfMat=z)
df_full3 = make_spikes(numUnits=100,rateProf=z[0,:])

plt.figure()
plt.plot(df_full3['timestamps'],df_full3['node_ids'],'.',markersize=1)

overallFR = 1000*df_full3.groupby('node_ids').count()/df_full3.timestamps.max()
plt.figure()
plt.subplot(2,1,1)
h = np.histogram(overallFR.values,bins=np.arange(0,16,1))
plt.bar(h[1][1:],h[0]/np.sum(h[0]))
plt.ylim(0,0.22)
plt.subplot(2,1,2)
g = np.histogram(np.exp(data.values),bins=np.arange(0,16,1))
plt.bar(g[1][1:],g[0]/np.sum(g[0]))
plt.ylim(0,0.22)

plt.figure()
plt.plot(g[1][1:],np.cumsum(g[0]/np.sum(g[0])),label='experimental')
plt.plot(h[1][1:],np.cumsum(h[0]/np.sum(h[0])),label='model')


ma_win = 50

#x1 = np.histogram(df_full['timestamps'],bins=np.arange(0,20000,1))
#y1 =  moving_average(x1[0],ma_win)

#x2 = np.histogram(df_full2['timestamps'],bins=np.arange(0,20000,1))
#y2 =  moving_average(x2[0],ma_win)

x3 = np.histogram(df_full3['timestamps'],bins=np.arange(0,20000,1))
y3 =  moving_average(x3[0],100)

plt.figure()

#plt.plot(x3[1][:-ma_win],moving_average(z[0,:],ma_win),label='input',color='k')
plt.plot(x3[1][:-100],y3*(1000/100),color='r',linestyle='-',label='output FR (100 units)')
plt.legend()
plt.ylabel('zscore firing rate')
plt.xlabel('time (ms)')

plt.show()

plt.figure()
plt.scatter(y*(1000/numUnits),moving_average(rateProf*1000,100),alpha=0.01)
plt.xlabel('actual FR (Hz)')
plt.ylabel('expected FR (Hz)')

#plt.plot((rateProf-np.mean(rateProf))/np.std(rateProf),label='input')
#plt.legend()

#plt.figure()
#plt.plot((y-np.mean(y))/np.std(y), (rateProf-np.mean(rateProf))/np.std(rateProf),'.')

##### CROSS CORRELATION #####
#def spkxcorr(spkref,spktgt,binwidth,window):
#    binwidth=1
#    window=100
#
#    h = np.zeros((1,window))
#    for i in spkref:
#        h += np.array(np.histogram(spktgt,bins=np.arange(i-int(window/2),i+int(window/2)+1,binwidth))[0])
#    return h
#
#spktgt = df_full.loc[df_full.node_ids==1,'timestamps']
#spkref = df_full.loc[df_full.node_ids==0,'timestamps']
#
#xcorrs = np.zeros((1000,100))
#for i in np.arange(1,1001):
#    xcorrs[i-1,:] = spkxcorr(spkref,spktgt,binwidth=1,window=100)
#
#
#
#
#plt.figure()
#plt.bar(np.arange(0,100),np.mean(xcorrs,0))



plt.show()
