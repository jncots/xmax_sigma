from abc import ABC, abstractmethod
import dataclasses


@dataclasses.dataclass
class CascadeParticle:
    pid: int
    energy: float
    xdepth: float

class ParticleEvent(ABC):
    
    def __init__(self, particle):
        self.set_particle(particle)

    @abstractmethod
    def set_particle(self, particle):
        pass
        
    @abstractmethod    
    def get_xdepth(self):
        pass
        
    @abstractmethod    
    def get_products(self):
        pass