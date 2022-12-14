import numpy as np
from MCEq import particlemanager


def rejection_sampler(p, xbounds, pmax):
    while True:
        x = np.random.rand(1) * (xbounds[1] - xbounds[0]) + xbounds[0]
        y = np.random.rand(1) * pmax
        if y <= p(x):
            return x


def rejection_sampler_rest(p, xbounds, pmax, xlimits):
    """Restricted rejection_sampler

    Args:
        p (_type_): _description_
        xbounds (_type_): _description_
        pmax (_type_): _description_
        xlimits (_type_): _description_

    Returns:
        _type_: _description_
    """
    while True:
        x = np.random.rand(1) * (xbounds[1] - xbounds[0]) + xbounds[0]
        if not (xlimits[0] <= x <= xlimits[1]):
            continue
        y = np.random.rand(1) * pmax
        if y <= p(x):
            return x


class MCEqDBFactory:
    """Class for creation of MCEqRun instance which is used
    as a database for particle decay information
    """

    from MCEq.core import MCEqRun
    import crflux.models as pm

    def get_db(self):
        return self.MCEqRun(
            interaction_model="SIBYLL2.3c",
            primary_model=(self.pm.HillasGaisser2012, "H3a"),
            theta_deg=0.0,
        )


class ParticleDecay:
    from particletools.tables import PYTHIAParticleData
    from scipy.interpolate import interp1d

    def __init__(self, mceq_db, pid, etot):
        self.mceq_db = mceq_db
        self.pythia_pdata = self.PYTHIAParticleData()
        self.pid = None
        self.etot = None
        self.set_decayed_particle(pid, etot)

    def set_decayed_particle(self, pid, etot):
        self.reset_pid_etot = False
        self.reset_mceq_data = False
        if self.pid != pid or self.etot != etot:
            self.reset_pid_etot = True
            self.reset_mceq_data = True

        if self.reset_pid_etot:
            self.pid = pid
            self.mass = self.pythia_pdata.mass(self.pid)
            self.etot = etot

            self.ekin = self.etot - self.mass
            gamma = self.etot / self.mass
            beta_gamma = np.sqrt((gamma + 1) * (gamma - 1))
            self.decay_length = beta_gamma * self.pythia_pdata.ctau(pid)

            self.reset_pid_etot = False

    def get_decay_length(self):
        """Returns decay length in `cm`
        for a particle set in constructor `ParticleDecay` or
        reset in `set_decayed_particle(pid, etot)` method
        """
        random_value = -np.log(1 - np.random.rand(1))
        return self.decay_length * random_value

    def get_decay_products(self):
        """Returns list of decay products `[(pid_1, etot_1),...]`
        for a particle set in constructor `ParticleDecay` or
        reset in `set_decayed_particle(pid, etot)` method
        """

        if self.reset_mceq_data:
            self._set_mceq_data()
            self.reset_mceq_data = False

        self.xleft = 1.0  # part of kinetic energy left for child particle

        # Get decay channel
        ind = np.argmax(np.random.multinomial(1, self.decay_probabilities, size=1))
        decay_products = self.decay_channels[ind][1]
        decay_size = len(decay_products)

        result = []

        for pid in decay_products:
            if decay_size > 1:
                result.append(self._get_child_particle(pid))
                decay_size -= 1
            else:
                result.append(self._get_last_child_particle(pid))

        return result

    def _set_mceq_data(self):
        self.parent_particle = particlemanager.MCEqParticle(
            self.pid, 0, self.mceq_db._energy_grid
        )
        self.parent_particle.set_decay_channels(self.mceq_db._decays, self.mceq_db.pman)
        # Set decay channels and probabilities
        self.decay_channels = self.pythia_pdata.decay_channels(self.pid)
        self.decay_probabilities = [channel[0] for channel in self.decay_channels]

    def _get_mceq_child_particle(self, pid):
        # !!Probably wrong!!
        # This is a workaround to get an accepted argument
        # for self.parent_particle.dNdec_dxlab
        # which can be different from what
        # self.pythia_pdata.decay_channels returns:
        # Example:
        # self.parent_particle.dNdec_dxlab accepts muon_left, muon_right
        # but self.pythia_pdata.decay_channels returns only muon
        # So this function returns first muon (muon_left or muon_right)
        # which it is found
        for child in self.parent_particle.decay_dists:
            if child.pdg_id[0] == pid:
                return self.mceq_db.pman[child.pdg_id]

    def _get_child_particle(self, child_pid):
        child_particle = self._get_mceq_child_particle(child_pid)
        data = self.parent_particle.dNdec_dxlab(self.ekin, child_particle)
        self.xgrid = data[0]
        self.pdf_max = np.max(data[1])
        self.pdf = self.interp1d(*data)

        # xpart is the part of kinetic energy of decaying particle,
        # i.e. E = xpart*ekin is the energy of child particle
        # we draw it from probability distribution function self.pdf
        xpart = rejection_sampler_rest(
            self.pdf,
            (self.xgrid[0], self.xgrid[-1]),
            self.pdf_max,
            (self.xgrid[0], self.xleft),
        )
 
        child_mass = self.pythia_pdata.mass(child_pid)
        # !!Probably wrong approach:
        # We use kinetic energy left from other particles
        # The correct way is to get it from joint pdf
        self.xleft -= xpart + child_mass / self.ekin
        return (child_pid, xpart * self.ekin + child_mass)

    def _get_last_child_particle(self, pid):
        return (pid, self.xleft * self.ekin + self.mass)


class DecayLength:
    from particletools.tables import PYTHIAParticleData

    def __init__(self, pdg, energy):
        pythia_pdata = self.PYTHIAParticleData()
        self.c_decay_time = pythia_pdata.ctau(pdg)
        self.mass = pythia_pdata.mass(pdg)
        self.set_energy(energy)

    def set_energy(self, energy):
        gamma = energy / self.mass
        beta_gamma = np.sqrt((gamma + 1) * (gamma - 1))
        self.decay_length = beta_gamma * self.c_decay_time

    def get_decay_length(self):
        random_value = -np.log(1 - np.random.rand(1))
        return self.decay_length * random_value
