import chromo
from particle_array import ParticleArray
import numpy as np


class HadronInteraction:
    def __init__(self):
        self.target =  (14, 7)
        ekin = chromo.kinematics.FixedTarget(20000, "proton", self.target)
        self.event_generator = chromo.models.Sibyll23d(ekin)
        # self.children = ParticleArray(100000)
        # self.failed_to_interact = None
    
    
    def run_event_generator(self, pstack):
        # pstack = ParticleArray()
        self.children = ParticleArray(100000)
        self.failed_to_interact = None
        
        for i in range(len(pstack)):
            try:
                self.event_generator.kinematics = chromo.kinematics.FixedTarget(
                pstack.energy[i], int(pstack.pid[i]), self.target
                )
            except ValueError:
                pstack.production_code[i] = 333
                
                    
            
            event = next(self.event_generator(1)).final_state()
            generation_num = pstack.generation_num[i] + 1
            self.children.push(pid = event.pid, 
                            energy = event.en, 
                            xdepth = np.full(len(event), pstack.xdepth_inter[i]),
                            generation_num = np.full(len(event), generation_num),
                            production_code = np.full(len(event), 1))
        
        
        
        self.failed_to_interact = pstack[np.where(pstack.production_code == 333)]
        pass    
        # print("small = ", pstack.pid[np.where(pstack.production_code == 333)])
        # print("small = ", pstack.energy[np.where(pstack.production_code == 333)])
        
        
    def get_children(self):
        return self.children

    def get_failed(self):
        return self.failed_to_interact

if __name__ == "__main__":
    
    hint = HadronInteraction()
    
    pstack = ParticleArray(1000)
    pstack.push(pid = 2212, energy = 1e5, xdepth = 0, xdepth_inter = 5, generation_num = 0)
    
    hint.run_event_generator(pstack)
    ch = hint.get_children().valid().copy()
    hint.run_event_generator(ch)
    ch = hint.get_children().valid().copy()
    hint.run_event_generator(ch)
    ch = hint.get_children().valid()
    print("pid = ", ch.pid)
    print("energy = ", ch.energy)
    print("xdepth = ", ch.xdepth)
    print("generation = ", ch.generation_num)
    print("valid_code = ", ch.valid_code)
    print(len(ch))
    