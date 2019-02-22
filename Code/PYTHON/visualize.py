import numpy as np
import matplotlib as mpl
from matplotlib import colors
import matplotlib.pyplot as plt

# define plot style

mpl.rcParams["font.family"]: monospace
mpl.rcParams["figure.subplot.wspace"]: 0.6
mpl.rcParams["figure.subplot.hspace"]: 0.6
mpl.rcParams["text.color"]: grey
mpl.rcParams["axes.titlesize"]: x-large
mpl.rcParams["axes.titlepad"]: 12.0
mpl.rcParams["axes.labelcolor"]: grey
mpl.rcParams["axes.labelpad"]: 8.0
mpl.rcParams["axes.spines.left"]: False
mpl.rcParams["axes.spines.top"]: False
mpl.rcParams["axes.spines.right"]: False
mpl.rcParams["axes.spines.bottom"]: False
mpl.rcParams["xtick.top"]: False
mpl.rcParams["xtick.color"]: grey
mpl.rcParams["ytick.right"]: False
mpl.rcParams["ytick.color"]: grey

# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))


def plot_ppc(data_dict, fit):
    # define variables
    Y = data_dict['Y'].flatten()
    Y_ppc = fit.extract(pars=['Y_ppc'])['Y_ppc']
    B = data_dict['A']+2
    
    fig = plt.figure(figsize=(8,8))
    ax = plt.gca()
    bins = [b - 0.5 for b in range(B + 1)]

    idxs = [ idx for idx in range(B) for r in range(2) ]
    xs = [ idx + delta for idx in range(B) for delta in [-0.5, 0.5]]

    # make a histogram for each sample of the markov iteration
    counts = [np.histogram(Y_ppc[n].flatten(), bins=bins)[0] for n in range(Y_ppc.shape[0])]
    probs = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    # find the percentiles for each sample-histogram
    creds = [np.percentile([count[b] for count in counts], probs) for b in range(B)]
    pad_creds = [ creds[idx] for idx in idxs ]

    ax.fill_between(xs, [c[0] for c in pad_creds], [c[8] for c in pad_creds], color='y', lw=0, alpha=.1)
    ax.fill_between(xs, [c[1] for c in pad_creds], [c[7] for c in pad_creds], color='orange', lw=0, alpha=.1)
    ax.fill_between(xs, [c[2] for c in pad_creds], [c[6] for c in pad_creds], color='orange', lw=0, alpha=.1)
    ax.fill_between(xs, [c[3] for c in pad_creds], [c[5] for c in pad_creds], color='red', lw=0, alpha=.1)

    ax.plot(xs, [c[4] for c in pad_creds], color='r', alpha=.3)
    ax.hist(Y, bins=bins, histtype='step', color='#4b0082', alpha=.8)
    plt.show()

    
def plot_betas(B1, B2):
    max_beta_01 = np.absolute(B1).max()
    max_beta_02 = np.absolute(B2).max()
    max_beta = min(max_beta_01, max_beta_02)
    # Plot the betas both generated and estimated
    plt.figure(figsize=(16,8))

    plt.subplot(211)
    plt.imshow(B1, cmap='RdGy_r', norm=MidpointNormalize(midpoint=0, vmin=-max_beta, vmax=max_beta))
    plt.colorbar()
    
    plt.subplot(212)
    plt.imshow(B2, cmap='RdGy_r', norm=MidpointNormalize(midpoint=0, vmin=-max_beta, vmax=max_beta))
    plt.colorbar()
    
    plt.show()

#     plt.subplot(312)
#     plt.plot(np.arange(12), B.mean(axis=1), color='r')
#     plt.title("Feature Level Avg")

#     plt.subplot(313)
#     y = B.mean(axis=0)
#     plt.plot(np.arange(len(y)), y, color='r')
#     plt.title("Respondent Avg")

def plot_fit_results(FIT, DATA, f=""):
    # plot results
    print("Results for fit {0}".format(f))
    B = FIT.extract(pars=['B'])['B'][-1].T
    utils.plot_betas(B, DATA['B'])
    utils.plot_ppc(DATA, FIT)


def plot_respondent(r, data_dict, fit):
    B = fit.extract(pars=['B'])['B']
    plt.figure(figsize=(12,15))
    for l in range(data_dict['L']):
        ax = plt.subplot(4,3,l+1)
        ax.hist(B[:,r,l], color='r', alpha=.5)
        ax.axvline(data_dict['B'][l,r], color='k')
        ax.set_title("Beta {0}".format(l+1))
    plt.show()

    plt.figure(figsize=(8,15))
    for t in range(data_dict['T']):
        ax = plt.subplot(5,2,t+1)
        ax.hist(Y_ppc[:,r,t], color='r', alpha=.5)
        ax.axvline(data_dict['Y'][r, t], color='k')
        ax.set_title("Task {0}".format(t+1))
    plt.show()


