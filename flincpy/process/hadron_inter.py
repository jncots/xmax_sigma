import chromo
from data_structs.particle_array import ParticleArray
import numpy as np


class HadronInteraction:
    def __init__(self):
        self.target = chromo.kinematics.CompositeTarget([("N", 0.78), ("O", 0.22)])
        # self.target = (14, 7)
        ekin = chromo.kinematics.FixedTarget(20000, "proton", self.target)
        self.event_generator = chromo.models.DpmjetIII191(ekin)
        # self.event_generator = chromo.models.EposLHC(ekin)
        # self.event_generator = chromo.models.Sibyll23d(ekin)
        # self.event_generator._ecm_min = 2 # GeV
        # chromo.models.Sibyll23d(ekin)
        # self.event_generator.set_unstable(111)
        # self.event_generator.set_unstable(-211)
        # self.event_generator.set_unstable(211)
    
    
    def run_event_generator(self, parents, children, failed_parents):
        
        number_of_interactions = 0
        pvalid = parents.valid()
        pvalid.production_code[:] = 0
        
        for i in range(len(pvalid)):
            try:
                self.event_generator.kinematics = chromo.kinematics.FixedTarget(
                    pvalid.energy[i], int(pvalid.pid[i]), self.target
                )
                
                if (self.event_generator.kinematics.ekin <= 2e0):
                    raise RuntimeError("Too low energy")
                
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
            # try:        
            #     event = next(self.event_generator(1)).final_state()
            # except RuntimeError as er:
            #     pvalid.production_code[i] = 333
            #     continue
                # print(er)
                # print(f"pdg_proj = {self.event_generator.kinematics.p1}" 
                #       f" and energy = {self.event_generator.kinematics.elab}")
                # print(self.event_generator.kinematics)   
                
            number_of_interactions += 1
            generation_num = pvalid.generation_num[i] + 1
            
            # if len(np.where(event.pid == 8)[0]) > 0:
            #     print(f"event.pid = {event.pid}")
            #     print(f"event.status = {event.status}")
            #     print(f"pdg_proj = {self.event_generator.kinematics.p1}" 
            #           f" and energy = {self.event_generator.kinematics.elab}")
            #     print(self.event_generator.kinematics) 
            #     self.event_generator._lib.dt_evtout(6)
            #     print("Here STOP")
            
            children.push(pid = event.pid, 
                            energy = event.en, 
                            xdepth = pvalid.xdepth_inter[i],
                            generation_num = generation_num,
                            parent_id = pvalid.id[i],
                            production_code = 777)
            
            # print("\n")
            # print(f"event.pid = {event.pid}")
            # print(f"event.pid = {event.en},\nsum = {np.sum(event.en)}, init = {pvalid.energy[i]}")
            # print(f"event.pid = {event.status}")
            # print(self.event_generator.kinematics)
            # print(f"Interaction: energy conservation {100*(np.sum(event.en) - pvalid.energy[i])/pvalid.energy[i]} %\n")
        
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
    
    print("failed pid = ", fp.pid)
    print("failed energy = ", fp.energy)
    print("failed xdepth = ", fp.xdepth)
    print("failed generation = ", fp.generation_num)
    print("failed valid_code = ", fp.production_code)
    print(len(fp))
    