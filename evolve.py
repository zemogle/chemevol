from scipy.integrate import quad
import functions as f
import numpy as np
from lookups import find_nearest
import logging
from lookups import t_lifetime

FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger('chem')


class ChemModel:
    def __init__(self, **inputs):
        #f.validate_initial_dict(inputs)
        try:
            self.gasmass_init = inputs['gasmass_init']
            self.gamma = inputs['gamma']
            self.imf_type = inputs['IMF_fn']
            self.dust_source = inputs['dust_source']
            self.destroy = inputs['destroy']
            self.inflows = inputs['inflows']
            self.outflows = inputs['outflows']
            self.SFH_file = inputs['SFH']
            if not self.SFH_file:
                self.SFH_file = 'Milkyway.sfh'
            self.sfh_file = self.SFH_file
            self.load_sfh()
        except KeyError:
            logger.error('You must provide initial parameters')

    def load_sfh(self):
        try:
            vals = np.loadtxt(self.SFH_file)
            self.sfh = vals
        except:
            logger.error("File '%s' will not parse" % self.SFH_file)
            self.sfh = None

    def sfr(self, t):
        try:
            vals = find_nearest(self.sfh,t)
            return vals[1]
        except:
            logger.error("No SFH yet")

    def ejected_mass(self, t, choice):
        mu = 120
        m = 0.8
        dm = 0.1
        em = 0.
        while m <= mu:
            m += dm
            em += f.ejected_gas_mass(m, self.sfr(t), choice) * dm
        return em

    def gas_mass(self, choice):
        t = self.sfh[0][0]
        t_end = self.sfh[-1][0]
        dlogt=(np.log10(t_end)-np.log10(t))/100.
        mg = 0.
        while t <= t_end:
            t += 10.**(np.log10(t)+dlogt)
            dmg = - self.sfr(t) + self.ejected_mass(t, choice) + f.inflows(self.sfr(t), self.inflows) + f.outflows(self.sfr(t), self.outflows)
            mg += dmg * 10.**dlogt
        return mg



def ejected_mass_integral(t):
    # integrate the ejected mass function

    return ej
