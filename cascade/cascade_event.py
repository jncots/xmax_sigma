from hadron_event import HadronEvent
from decay_event import DecayEvent
from xdepth_conversion import XdepthConversion




class CascadeEvent:
    def __init__(self, emin_threshold, particle):
        
        self.emin_threshold = emin_threshold
        self.particle = particle
        self.hadron_event = HadronEvent(self.particle)
        self.decay_event = DecayEvent(self.particle)
        
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
    # def _get_event(self, particle):
    #     """Without decay. Only for debugging

    #     Args:
    #         particle (_type_): _description_

    #     Returns:
    #         _type_: _description_
    #     """
    #     self.hadron_event.set_particle(particle)
    #     xdepth_hadron = self.hadron_event.get_xdepth()
        
    #     if not xdepth_hadron:
    #         return None
    #     return self.hadron_event.get_products()
    
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
                
    def get_event_particles(self, cur_particle):
        
        event = self._get_event(cur_particle)
        
        if not event:
            return None
        
        inter_particles = []
        final_particles = []
        
        for particle in event:
            if (particle.energy < self.emin_threshold) or (particle.pid not in self.hadron_event.valid_pids):
                final_particles.append(particle)
            else:
                inter_particles.append(particle)
        return inter_particles, final_particles