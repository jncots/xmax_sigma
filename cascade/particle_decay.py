import numpy as np
import random
from pythia_decay import Pythia8DecayAfterburner
from particletools.tables import PYTHIAParticleData
from xdepth_conversion import XdepthConversion


class ParticleDecay:
    def __init__(self):
        """Initialization

        Args:
            pid (int): pdg id of particle
            etot (float): total energy of particle in GeV
        """

        self.decay_event = Pythia8DecayAfterburner()
        self.pythia_pdata = PYTHIAParticleData()
        self.xconv = XdepthConversion()
        self.max_height = self.xconv.get_max_height()

    def get_decay_length(self, particle):
        """Returns decay length in `cm`"""

        ctau0 = float(self.pythia_pdata.ctau(particle.pid))
        if ctau0 == 0 or ctau0 == np.inf:
            return None

        gamma = particle.energy / float(self.pythia_pdata.mass(particle.pid))
        # mean_length = np.sqrt((gamma + 1) * (gamma - 1)) * ctau0
        # random_value = -np.log(1 - random.random())
        # return random_value * mean_length
        return -np.log(1 - random.random()) * np.sqrt((gamma + 1) * (gamma - 1)) * ctau0

    def get_xdepth(self, particle):
        length = self.get_decay_length(particle)
        if length is None:
            return None
        return self.xconv.get_delta_xdepth(particle.xdepth, length)
    
    def set_xdepth_decay(self, particle):
        xdepth = self.get_xdepth(particle)
        
        if xdepth is None:
            particle.xdepth_decay = None
        else:    
            particle.xdepth_decay = particle.xdepth + xdepth

    def get_decay_products(self, particle_list):
        return self.decay_event(particle_list)


if __name__ == "__main__":
    from particle_event import CascadeParticle
    particle = CascadeParticle(111, 1e10, 100)
    dp = ParticleDecay()
    print(dp.get_xdepth(particle))
    