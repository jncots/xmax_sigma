from hadron_event import HadronEvent
from decay_event import DecayEvent




class CascadeEvent:
    def __init__(self, emin_threshold, particle):
        
        self.emin_threshold = emin_threshold
        self.particle = particle
        self.hadron_event = HadronEvent(self.particle)
        self.decay_event = DecayEvent(self.particle)

    def _get_event(self, particle):
        
        self.hadron_event.set_particle(particle)
        self.decay_event.set_particle(particle)
        
        xdepth_hadron = self.hadron_event.get_xdepth()
        xdepth_decay = self.decay_event.get_xdepth()

        if xdepth_decay and xdepth_decay < xdepth_hadron:
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