import numpy as np
import scipy.constants as sc

# Each entry: 'param_name': ('prior_type', [prior_params])

# Prior types:
#   'normal'    : [mean, std]        — Gaussian prior
#   'uniform'   : [low, high]        — Flat prior between bounds
#   'linewidth' : ???????            — Physical lower bound from thermal
#                                      width + upper bound you set [m/s]

priors_dict = {
    'inc'       : ('normal',    [30.7, 0.7]),   # inclination (deg)
    'PA'        : ('normal',    [109.0,1.8]),   # position angle (deg)
    'mstar'     : ('uniform',   [0.44, 1.5]),   # stellar mass (Msun)
    'Rout'       : ('normal',    [60, 30]),   # outer radius (au)
    'Tb_10'     : ('uniform',   [20, 200]),   # brightness temp at 10 au (K)
    'Tb_q'      : ('normal',    [-0.7, 0.1]),   # brightness temp gradient
    'tau_0'    : ('uniform',   [1.0, 30.0]),   # log10 optical depth at 10 au
    'vlsr'      : ('normal',    [1.2e4, 5e2]),   # systemic velocity (m/s)
    'dx'        : ('normal',    [0.0, 0.1]),   # RA offset (arcsec)
    'dy'        : ('normal',    [0.0, 0.1]),   # Dec offset (arcsec)
}

# --- Unpack into lists that csalt expects (do not edit below this line) ---
pri_types = [v[0] for v in priors_dict.values()]
pri_pars  = [v[1] for v in priors_dict.values()]

# --- Pre-defined prior functions ---

def logprior_uniform(theta, ppars):
    """Flat prior between ppars[0] and ppars[1]."""
    if np.logical_and((theta >= ppars[0]), (theta <= ppars[1])):
        return np.log(1 / (ppars[1] - ppars[0]))
    else:
        return -np.inf

def logprior_normal(theta, ppars):
    """Gaussian prior with mean ppars[0] and std ppars[1]."""
    return (-np.log(ppars[1] * np.sqrt(2 * np.pi))
            - 0.5 * ((theta - ppars[0])**2 / ppars[1]**2))

def logprior_linewidth(theta, ppars):
    lw0 = np.sqrt(2 * sc.k * ppars[1] / (28 * (sc.m_p + sc.m_e)))
    if np.logical_and((theta >= lw0), (theta <= ppars[0])):
        return 0
    else:
        return -np.inf

# --- Log-prior calculator  ---

def logprior(theta):
    logptheta = np.empty_like(theta)
    # linewidth prior needs Tmax_b (index 8) as second parameter
    #pri_pars[9] = [pri_pars[9][0], theta[8]]
    for i in range(len(theta)):
        cmd = ('logprior_' + pri_types[i] +
               '(theta[' + str(i) + '], ' + str(pri_pars[i]) + ')')
        logptheta[i] = eval(cmd)
    return logptheta