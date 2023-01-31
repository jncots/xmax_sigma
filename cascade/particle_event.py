from abc import ABC, abstractmethod
import dataclasses
from typing import List, Tuple


@dataclasses.dataclass
class CascadeParticle:
    pid: int
    energy: float
    xdepth: float
    production_mode: int = 0 
    # 0 default, 
    # 1 interaction, 
    # 2 decay
    generation_number: int = 0
    final_code: int = 0 
    # 0 default - may interact (not final)
    # 1 interaction not supported by generator, 
    # 2 below threshold
    # 3 interaction point below ground (xdepth > max xdepth)
    xdepth_decay: float = 0
    parent: List = dataclasses.field(default_factory=list)
    
    def get_parents(self):
        parents = []
        
        parent = self.parent[0] if self.parent else None
        while parent:
            parents.append((parent.pid, parent.energy, parent.xdepth, parent.production_mode))
            parent = parent.parent[0] if parent.parent else None
        
        ngen = len(parents) - 1
        results = {}
        for i, p in enumerate(parents):
            results[ngen - i] = p
        return results        
        


class ParticleEvent(ABC):
    def __init__(self, particle):
        self.max_xdepth = 10000
        self.set_particle(particle)

    def set_max_xdepth(self, max_xdepth):
        self.max_xdepth = max_xdepth

    @abstractmethod
    def set_particle(self, particle):
        pass

    @abstractmethod
    def get_xdepth(self):
        pass

    @abstractmethod
    def get_products(self):
        pass
