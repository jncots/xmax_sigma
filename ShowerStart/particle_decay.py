import numpy as np
import crflux.models as pm
from MCEq import particlemanager



def rejection_sampler(p,xbounds,pmax):
    while True:
        x = np.random.rand(1)*(xbounds[1]-xbounds[0])+xbounds[0]
        y = np.random.rand(1)*pmax
        if y<=p(x):
            return x

def restricted_rej_sampler(p, xbounds, pmax, xlimits):
    while True:
        x = np.random.rand(1)*(xbounds[1]-xbounds[0])+xbounds[0]
        if not (xlimits[0] <= x <= xlimits[1]):
            continue
        y = np.random.rand(1)*pmax
        if y<=p(x):
            return x

class ParticleDecay:
    from MCEq.core import MCEqRun
    from particletools.tables import PYTHIAParticleData
    from scipy.interpolate import interp1d
    
    def __init__(self):
        self.pythia_pdata = self.PYTHIAParticleData()
        self.mceq_run = self.MCEqRun(
            interaction_model='SIBYLL2.3c',
            primary_model = (pm.HillasGaisser2012, "H3a"),
            theta_deg=0.0)
        
    def set_decayed_particle(self, pid):
        mcr = self.mceq_run    
        self.particle = particlemanager.MCEqParticle(pid, 0, mcr._energy_grid)
        self.particle.set_decay_channels(mcr._decays, mcr.pman)
        self.particle_mass = self.pythia_pdata.mass(pid)
        
        # Set decay channels
        self.decay_channels = self.pythia_pdata.decay_channels(pid)
        self.decay_probabilities = [channel[0] for channel in self.decay_channels] 
    
    def _to_ekin(self, etot):
        return etot - self.particle_mass
    
    def _to_etot(self, ekin):
        return ekin + self.particle_mass
    
    def _child_particle_pdf(self, etot, pid):
        ekin = self._to_ekin(etot)
        data = self.particle.dNdec_dxlab(ekin, self.mceq_run.pman[pid])
        self.particle_ekin = ekin
        self.xgrid = data[0]
        self.pdf_max = np.max(data[1])
        self.pdf = self.interp1d(*data)
        
    def get_decay_channel(self):
        ind = np.argmax(np.random.multinomial(1, self.decay_probabilities, size=1))
        return self.decay_channels[ind]
    
    
    def get_child_particle_energy(self, xmax):
        self._child_particle_pdf()
        
        xpart = restricted_rej_sampler(self.pdf,
                                  (self.xgrid[0], self.xgrid[-1]),
                                  self.pdf_max,
                                  (self.xgrid[0], xmax))
        return xpart