from hadron_event import HadronEvent
from decay_event import DecayEvent
from xdepth_conversion import XdepthConversion
from pythia_decay import DecayByPythia
from particle_event import CascadeParticle




class CascadeEvent:
    def __init__(self, emin_threshold, particle):
        
        self.emin_threshold = emin_threshold
        self.particle = particle
        self.hadron_event = HadronEvent(self.particle)
        self.decay_event = DecayEvent(self.particle)
        self.pythia_dec = DecayByPythia()
        
        max_xdepth = self._get_max_depth()
        print(f"CascadeEvent max_xdepth = {max_xdepth}")     
        self.hadron_event.set_max_xdepth(max_xdepth)
        self.decay_event.set_max_xdepth(max_xdepth)

        self._decay_on = True

    def _get_max_depth(self):
        xconv = XdepthConversion()
        return xconv.get_max_xdepth()
        
    
    def _reset_ncalls(self):
        self.hadron_event._get_prod_ncalls = 0
        self.decay_event._get_prod_ncalls = 0
    
    def set_decay_on(self, decay_on):
        self._decay_on = decay_on

    def _get_event(self, particle):
        
        
        if self._decay_on:
            self.decay_event.set_particle(particle)
            xdepth_decay = self.decay_event.get_xdepth()
        else:
            xdepth_decay = None    
               
        self.hadron_event.set_particle(particle)
        xdepth_hadron = self.hadron_event.get_xdepth()
        
        if (not xdepth_decay) and (not xdepth_hadron):
            return None

        if xdepth_decay and (not xdepth_hadron):
            return self.decay_event.get_products()
        
        if (not xdepth_decay) and xdepth_hadron:
            return self.hadron_event.get_products()
        
        if xdepth_decay and xdepth_hadron:
            if xdepth_decay < xdepth_hadron:
                return self.decay_event.get_products()
            else:
                return self.hadron_event.get_products()
            
    def decay_particles(self, prod):
        
        result = []
        for p in prod:
            
            dprod = self.pythia_dec.get_decayed_products(p.pid, p.energy)
            if dprod:
                for pp in dprod:
                    cp = CascadeParticle(pp[0], pp[1], p.xdepth, 2, p.generation_number + 1)
                    result.append(cp)
            else:
                result.append(p)    
        return result
                         
                
    def get_event_particles(self, cur_particle):
                
        event = self._get_event(cur_particle)
        
        if not event:
            return None
        
        inter_particles = []
        final_particles = []
        
        for particle in event:
            if (particle.energy < self.emin_threshold):
                final_particles.append(particle)
            else:
                inter_particles.append(particle)
                
        final_particles = self.decay_particles(final_particles)        
                
        return inter_particles, final_particles