import numpy as np
import scipy.constants as sc

# Each entry: 'param_name': ('prior_type', [prior_params])

# Prior types:
#   'normal'    : [mean, std]        — Gaussian prior
#   'uniform'   : [low, high]        — Flat prior between bounds
#   'linewidth' : ???????            — Physical lower bound from thermal
#                                      width + upper bound you set [m/s]

priors_dict = {
    'inc'       : ('normal',    [65.8, 0.2 ]),   # inclination (deg)
    'PA'        : ('normal',    [324.9,0.3]),   # position angle (deg)
    'mstar'     : ('uniform',   [0.2, 1.5]),   # stellar mass (Msun)
    'r_l'       : ('uniform',   [150, 50]),   # outer radius (au)
    'z_10'      : ('uniform',   [0.0, 5.0]),   # emission height norm (" at r=1")
    'z_q'       : ('normal',    [1.25, 0.5]),   # emission height gradient
    'Tb_10'     : ('uniform',   [20, 200]),   # brightness temp at 10 au (K)
    'Tb_q'      : ('normal',    [-0.5, 0.2]),   # brightness temp gradient
    'Tmax_b'    : ('uniform',   [5, 100]),   # max brightness temp (K)
    'dV_10'     : ('linewidth', [1000]),    # linewidth at 10 au (m/s)
    'logtau_10' : ('uniform',   [0.0, 2.0]),   # log10 optical depth at 10 au
    'tau_q'     : ('uniform',   [-4, 1]),   # optical depth gradient
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
    pri_pars[9] = [pri_pars[9][0], theta[8]]
    for i in range(len(theta)):
        cmd = ('logprior_' + pri_types[i] +
               '(theta[' + str(i) + '], ' + str(pri_pars[i]) + ')')
        logptheta[i] = eval(cmd)
    return logptheta