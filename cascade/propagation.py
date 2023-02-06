from hadron_event import HadronEvent
from particle_decay import ParticleDecay
from particle_array import ParticleArray

class Propagation:
    
    def __init__(self):
        self.hadron_event = HadronEvent()
        self.decay_event = ParticleDecay()
        
    
    def propagate(self, pstack):            
        self.calc_xd_decay(pstack.dxdepth_decay[0:len(pstack)]) 
        
        # pstack.dxdepth_decay[0:len(pstack)] = -np.log(1 - np.random.rand(len(pstack)))
        self.calx_xd_inter(pstack.dxdepth_inter[0:len(pstack)])   
    
    
    def calc_xd_decay(self, xdepth_decay):
        
        xdepth_decay[:] = -np.log(1 - np.random.rand(len(xdepth_decay)))  
        
        
        
    def calx_xd_inter(self, xdepth_inter):
        xdepth_inter[:] = -np.log(1 - np.random.rand(len(xdepth_inter)))
        # print(self.decay_event.get_xdepth())
        
        
if __name__ == "__main__":
    import numpy as np
    pstack = ParticleArray(10)
    
    pstack.push(pid = np.array([11, 12, 13, 14]), energy = np.array([100, 101, 102, 103]))
    print(pstack.pid[0:len(pstack)])
    print(pstack.energy[0:len(pstack)])
    print(pstack.dxdepth_decay[0:len(pstack)])
    print(pstack.dxdepth_inter[0:len(pstack)])
    
    Propagation().propagate(pstack)
    
    print("After propagation")
    print(pstack.dxdepth_decay[0:len(pstack)])
    print(pstack.dxdepth_inter[0:len(pstack)])