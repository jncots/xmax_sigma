
from particle_array import ParticleArray

class CascadeDriver:
    def __init__(self):
        self.wstack = ParticleArray(10000)
        self.fstack = ParticleArray(10000)
    
    
    def run(self, pdg, energy, xdepth = 0):
        self.wstack.push(pid = pdg, 
                         energy = energy, 
                         xdepth = xdepth)
        
        self.fstack 
        
    
    
    
    def get_final_particles(self):
        return self.fstack
        
    
    
    
    