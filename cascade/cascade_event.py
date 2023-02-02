from hadron_event import HadronEvent
from particle_decay import ParticleDecay
from pythia_decay import Pythia8DecayAfterburner
from particle_event import CascadeParticle




class CascadeEvent:
    def __init__(self, emin_threshold):
        
        self.emin_threshold = emin_threshold
        self.hadron_event = HadronEvent()
        self.decay_event = ParticleDecay()
        self.afterburner = Pythia8DecayAfterburner()
        self.set_decay_on(True)
        
        self.interacting_particles = []
        self.decaying_particles = []
        self.final_particles = []
        
        self.stable_particle_pids = [2212, 2112, 22, 11, 12, 13, 14]
        
    
    def _reset_ncalls(self):
        self.hadron_event._get_prod_ncalls = 0
        self.decay_event._get_prod_ncalls = 0
    
    def set_decay_on(self, decay_on = True):
        self._decay_on = decay_on
            
    def decay_particles(self, prod):
        return self.afterburner(prod)
        
    def run_event(self, particle):
        
        xdepth_decay = self.decay_event.get_xdepth(particle)
        xdepth_hadron = self.hadron_event.get_xdepth(particle)
        
        if xdepth_decay is None: 
            if xdepth_hadron is None:
                self.run_empty_event(particle)
            else:
                self.run_hadron_event(particle, xdepth_hadron)
        else:
            if xdepth_hadron is None:
                self.run_decay_event(particle, xdepth_decay)
            else:
                if xdepth_decay < xdepth_hadron:
                    self.run_decay_event(particle, xdepth_decay)
                else:
                    self.run_hadron_event(particle, xdepth_hadron)
                   
                    
        
    
    def _clear_buffers(self):    
        self.interacting_particles = []
        self.decaying_particles = []
        self.final_particles = []
        
        
    def run_empty_event(self, particle):
        self._clear_buffers()
        particle.final_code = 1
        self.final_particles.append(particle)
        
                
    def run_hadron_event(self, particle, xdepth_hadron):
        self._clear_buffers()
        
        products = self.hadron_event.get_products(particle, xdepth_hadron)
        if products is None:
            self.run_empty_event(particle)
            return
        
        for particle in products:            
            if (particle.energy < self.emin_threshold):
                if abs(particle.pid) not in self.stable_particle_pids:
                    self.decay_event.set_xdepth_decay(particle)
                    
                    if particle.xdepth_decay is None:
                        particle.final_code = 3
                        particle.xdepth = self.hadron_event.xdepth_max
                        self.final_particles.append(particle)
                    else:    
                        self.decaying_particles.append(particle)
                else:       
                    particle.final_code = 2
                    self.final_particles.append(particle)
            else:
                self.interacting_particles.append(particle)
                    
    
    def run_decay_event(self, particle, xdepth_decay):
        # print(f"Decaying particle: {particle.pid}")
        self._clear_buffers()
        particle.xdepth_decay = particle.xdepth + xdepth_decay
        self.decaying_particles.append(particle)
        
                             
                
    def get_event_particles(self, particle):   
        self.run_event(particle)
        return self.interacting_particles, self.final_particles, self.decaying_particles
    
    

if __name__ == "__main__":
    cev = CascadeEvent(1e4)
    particle = CascadeParticle(2212, 1e8, 0)
    
    
    res = cev.get_event_particles(particle)
    print(len(res[0]), len(res[1]), len(res[2]))
    
    
        