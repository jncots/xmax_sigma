import chromo
from particle_array import ParticleArray
import numpy as np
import logging


class HadronInteraction:
    def __init__(self):
        self.target = (14, 7)
        ekin = chromo.kinematics.FixedTarget(20000, "proton", self.target)
        self.event_generator = chromo.models.Sibyll23d(ekin)
    
    
    def run_event_generator(self, parents, children, failed_parents):
        
        number_of_interactions = 0
        pvalid = parents.valid()
        pvalid.production_code[:] = 0
        
        for i in range(len(pvalid)):
            try:
                self.event_generator.kinematics = chromo.kinematics.FixedTarget(
                    pvalid.energy[i], int(pvalid.pid[i]), self.target
                )
            except Exception as e:
                if "projectile" in str(e):
                    # projectile is not allowed
                    pvalid.production_code[i] = 111
                elif "center-of-mass" in str(e):
                    # center-of-mass energy  < minimum energy 10.0 GeV
                    pvalid.production_code[i] = 222
                else:
                    pvalid.production_code[i] = 333
                
                continue    
                    
            event = next(self.event_generator(1)).final_state()
            number_of_interactions += 1
            generation_num = pvalid.generation_num[i] + 1
            children.push(pid = event.pid, 
                            energy = event.en, 
                            xdepth = pvalid.xdepth_inter[i],
                            generation_num = generation_num,
                            production_code = 777)
        
        
        
        failed_parents.append(pvalid[np.where(pvalid.production_code > 0)])
        
        return number_of_interactions
        

if __name__ == "__main__":
    
    hint = HadronInteraction()
    
    parents = ParticleArray()
    children = ParticleArray()
    failed_parents = ParticleArray()
    
    parents.push(pid = 2212, energy = 1e5, xdepth = 0, xdepth_inter = 5, generation_num = 0)
    
    hint.run_event_generator(parents, children, failed_parents)
    parents1 = children.valid().copy()
    hint.run_event_generator(parents1, children, failed_parents)
    # hint.run_event_generator(ch)
    # ch = hint.get_children().valid().copy()
    # hint.run_event_generator(ch)
    # ch = hint.get_children().valid()
    ch = children.valid()
    fp = failed_parents.valid()
    print("pid = ", ch.pid)
    print("energy = ", ch.energy)
    print("xdepth = ", ch.xdepth)
    print("generation = ", ch.generation_num)
    print("valid_code = ", ch.valid_code)
    print(len(ch))
    
    print("pid = ", fp.pid)
    print("energy = ", fp.energy)
    print("xdepth = ", fp.xdepth)
    print("generation = ", fp.generation_num)
    print("valid_code = ", fp.production_code)
    print(len(fp))
    