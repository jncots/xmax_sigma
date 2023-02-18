from decay_xdepth import DecayXdepth
from csec_tables import CrossSectionTableMCEq
from pdg_pid_map import PdgLists
from csec_on_table import CrossSectionOnTable
from xdepth_on_table import XdepthOnTable

class NextDecayXdepth:
    def __init__(self, xdepth_on_table):
        self.decay_xdepth = DecayXdepth(xdepth_on_table)
         
    def get_xdepth(self, pstack):        
        pslice = slice(0, len(pstack))
        pdg = pstack.pid[pslice]
        energy = pstack.energy[pslice]
        xdepth = pstack.xdepth[pslice]
        pstack.xdepth_decay[pslice] = (self.decay_xdepth
                                       .get_xdepth(pdg, energy, xdepth))
        
class NextInterXdepth:
    def __init__(self, xdepth_on_table):
        cs_table = CrossSectionTableMCEq()
        cs_table.add_pdgs(PdgLists().longer_pi0_to_mceq)
        self.inter_xdepth = CrossSectionOnTable(cs_table)
        self.xdepth_max = xdepth_on_table
         
    def get_xdepth(self, pstack):       
        pslice = slice(0, len(pstack))
        pdg = pstack.pid[pslice]
        energy = pstack.energy[pslice]
        xdepth = pstack.xdepth[pslice]
        pstack.xdepth_inter[pslice] = (self.inter_xdepth
                                       .get_xdepth(pdg, energy) + xdepth)       
        
        
        
if __name__ == "__main__":
    from particle_array import ParticleArray
    import numpy as np
    
    next_decay = NextDecayXdepth()
    next_inter = NextInterXdepth()
    
    pstack = ParticleArray(size=100)
    pstack.push(
        pid=np.array([111, 22, 111, 13, -13, 2212]),
        energy=np.array([2e10, 2e0, 2e10, 1e0, 1e0, 1e0]),
        xdepth=np.array([100, 56, 100, 98, 56, 500]),
    )
    
    print(pstack.xdepth_decay[0:len(pstack)])
    print(pstack.xdepth_inter[0:len(pstack)])
    
    next_decay.get_xdepth(pstack)
    next_inter.get_xdepth(pstack)
    
    print("After")
    print(pstack.xdepth_decay[0:len(pstack)])
    print(pstack.xdepth_inter[0:len(pstack)])
    
    
         