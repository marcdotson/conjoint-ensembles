import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.special import erfinv
import matplotlib.pyplot as plt

C = 2
K = 2
J = 100
S = 5
G = 1

# impulse priors
alpha = 1

T = np.random.uniform(0,1,size=10000)

draws_from_l1_prior = np.array([-1/alpha * np.log(1-t) for t in T])
draws_from_cauchy_prior = np.array([1/alpha * np.tan(np.pi/2 * t) for t in T])
draws_from_white_noise_prior = np.array([alpha*np.sqrt(2)*erfinv(t) for t in T])

plt.figure(figsize=(12,8))
plt.subplot(131)
plt.imshow(draws_from_l1_prior.reshape(100,100),cmap='Greys_r')
plt.subplot(132)
plt.imshow(draws_from_cauchy_prior.reshape(100,100),cmap='Greys_r')
plt.subplot(133)
plt.imshow(draws_from_white_noise_prior.reshape(100,100),cmap='Greys_r')
plt.show()
