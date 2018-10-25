############################################################
# Initial setup
############################################################

import matplotlib.pyplot as plt
import numpy as np
import pystan
import stan_utility

############################################################
# Fit Poisson model
############################################################
data = pystan.read_rdump('discrete_regression.data.R')
data['theta'] = 10
model = stan_utility.compile_model('poisson2.stan')
fit = model.sampling(data=data) #, seed=4938483)
params = fit.extract()

# Perform a posterior predictive check by plotting
# the posterior predictive distribution of various
# components against the observed data

plt.figure(figsize=(16,8))
ax = plt.subplot(121)
B = 21

bins = [b - 0.5 for b in range(B + 1)]

idxs = [ idx for idx in range(B) for r in range(2) ]
xs = [ idx + delta for idx in range(B) for delta in [-0.5, 0.5]]

# make a histogram for each sample of the markov iteration
counts = [np.histogram(params['y_ppc'][n], bins=bins)[0] for n in range(4000)]
probs = [10, 20, 30, 40, 50, 60, 70, 80, 90]
# find the percentiles for each sample-histogram
creds = [np.percentile([count[b] for count in counts], probs) for b in range(B)]
pad_creds = [ creds[idx] for idx in idxs ]

ax.fill_between(xs, [c[0] for c in pad_creds], [c[8] for c in pad_creds], color='grey', lw=0, alpha=.1)
ax.fill_between(xs, [c[1] for c in pad_creds], [c[7] for c in pad_creds], color='grey', lw=0, alpha=.1)
ax.fill_between(xs, [c[2] for c in pad_creds], [c[6] for c in pad_creds], color='grey', lw=0, alpha=.1)
ax.fill_between(xs, [c[3] for c in pad_creds], [c[5] for c in pad_creds], color='grey', lw=0, alpha=.1)

ax.plot(xs, [c[4] for c in pad_creds], color='grey', alpha=.2)
ax.hist(data['y'], bins=bins, histtype='step', lw=3, color='#4b0082', alpha=.8)

ax.set_xlabel("y")
ax.set_ylabel("Posterior Predictive Distribution")
ax.set_title("Poisson")
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_ylim((0,70))

### MY TURN ###

model = stan_utility.compile_model('poisson1.stan')
fit = model.sampling(data=data) #, seed=4938483)
params = fit.extract()

ax = plt.subplot(122)
B = 21

bins = [b - 0.5 for b in range(B + 1)]

idxs = [ idx for idx in range(B) for r in range(2) ]
xs = [ idx + delta for idx in range(B) for delta in [-0.5, 0.5]]

# make a histogram for each sample of the markov iteration
counts = [np.histogram(params['y_ppc'][n], bins=bins)[0] for n in range(4000)]
probs = [10, 20, 30, 40, 50, 60, 70, 80, 90]
# find the percentiles for each sample-histogram
creds = [np.percentile([count[b] for count in counts], probs) for b in range(B)]
pad_creds = [ creds[idx] for idx in idxs ]

ax.fill_between(xs, [c[0] for c in pad_creds], [c[8] for c in pad_creds], color='y', lw=0, alpha=.1)
ax.fill_between(xs, [c[1] for c in pad_creds], [c[7] for c in pad_creds], color='orange', lw=0, alpha=.1)
ax.fill_between(xs, [c[2] for c in pad_creds], [c[6] for c in pad_creds], color='orange', lw=0, alpha=.1)
ax.fill_between(xs, [c[3] for c in pad_creds], [c[5] for c in pad_creds], color='red', lw=0, alpha=.1)

ax.plot(xs, [c[4] for c in pad_creds], color='r', alpha=.3)
ax.hist(data['y'], bins=bins, histtype='step', color='#4b0082', alpha=.8)

ax.set_xlabel("y")
ax.set_ylabel("Posterior Predictive Distribution")
ax.set_title("Negative Binomial")
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_ylim((0,70))


plt.show()
