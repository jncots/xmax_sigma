import numpy as np
from particle_array import ParticleArray, FilterCode
from pdg_pid_map import PdgLists
from particle_xdepths import default_xdepth_getter
from hadron_inter import HadronInteraction
from decay_driver import DecayDriver
import numpy as np

class CascadeDriver:
    def __init__(self, xdepth_getter = default_xdepth_getter, threshold_energy = 1e3):
        self.threshold_energy = threshold_energy
        self.xdepth_getter = xdepth_getter
        self.pdg_lists = PdgLists()
        self.hadron_interaction = HadronInteraction()
        self.decay_driver = DecayDriver(default_xdepth_getter, 
                            decaying_pdgs=[111], 
                            stable_pdgs=[-211, 211, -13, 13],
                            # decaying_pdgs=[111, -211, 211, -13, 13],
                            )
        
        
        self.working_stack = ParticleArray()
        self.final_stack = ParticleArray()
        self.decay_stack = ParticleArray()
        self.inter_stack = ParticleArray()
        self.failed_stack = ParticleArray()
        self.children_stack = ParticleArray()
    
    
    def run(self, *, pdg, energy, xdepth = 0):
        
        self.number_of_decays = 0
        self.number_of_interactions = 0
        
        self.working_stack.push(pid = pdg, 
                         energy = energy, 
                         xdepth = xdepth)
        
        iloop = 1
        while len(self.working_stack) > 0:
            
            print(f"{iloop} Number of inter = {self.number_of_interactions}"
                  f" number of decays = {self.number_of_decays}")
            self.set_xdepths()
            self.run_interaction()
            self.filter_particles()
            
            if len(self.working_stack) == 0:
                self.decay_particles()
            
            iloop += 1    
        #     self.cascade_event.get_event_particles(self.working_stack)
        
        # self.final_stack.append(self.working_stack)
        
    # def filter_below_threshold(self):
    #     threshold_energy = 1e3
        
    #     wstack = self.working_stack.valid()
    #     finals_slice = np.where(wstack.energy < threshold_energy & wstack.valid_code == 1)[0]
        
    #     self.decay_stack.append(wstack[finals_slice])
    #     wstack.valid_code[finals_slice] = 0
        
    
    
    def set_xdepths(self):
        """pstack is assumed to contain particles which:
            * has energy > threshold_energy
            * xdepth < surface_xdepth
            * has pdg which can be can decay or interact

        Args:
            pstack (_type_): _description_

        Returns:
            _type_: _description_
        """
        
        self.xdepth_getter.get_decay_xdepth(self.working_stack)
        self.xdepth_getter.get_inter_xdepth(self.working_stack)
        
        
        pvalid = self.working_stack.valid()
        max_xdepth = self.xdepth_getter.max_xdepth
        
        # Copies of particles with specific properties
        self.final_stack.append(pvalid[np.where((pvalid.xdepth_inter >= max_xdepth) &
                                (pvalid.xdepth_decay >= max_xdepth))[0]])
        
        self.inter_stack = pvalid[np.where(pvalid.xdepth_inter < pvalid.xdepth_decay)[0]]
        self.decay_stack.append(pvalid[np.where(pvalid.xdepth_decay < pvalid.xdepth_inter )[0]])
        
    
    def run_interaction(self):
        self.children_stack.clear()
        self.failed_stack.clear()
        self.number_of_interactions += self.hadron_interaction.run_event_generator(
                                                    parents = self.inter_stack, 
                                                    children = self.children_stack, 
                                                    failed_parents = self.failed_stack)
              
        if len(self.failed_stack) > 0:
            fvalid = self.failed_stack.valid()
            finals_true = np.isin(fvalid.pid, self.pdg_lists.mceq_finals)
            finals_false = np.logical_not(finals_true)           
            self.final_stack.append(fvalid[np.where(finals_true)])
            decaying_particles = fvalid[np.where(finals_false)]
            decaying_particles.filter_code[:] = FilterCode.XD_DECAY_OFF.value
            self.decay_stack.append(decaying_particles)
           
         
    def filter_particles(self):
        decay_pdgs = np.array([111], dtype = np.int32)
        
        cstack = self.children_stack.valid()
        
        decaying_true = np.isin(cstack.pid, decay_pdgs)
        self.decay_stack.append(cstack[np.where(decaying_true)[0]])
        
        finals_true = np.logical_or(np.isin(cstack.pid, self.pdg_lists.mceq_finals), 
                                    (cstack.energy < self.threshold_energy))
        self.final_stack.append(cstack[np.where(finals_true)[0]])
        
        working_true = np.logical_not(np.logical_or(decaying_true, finals_true))
        self.working_stack.clear()
        self.working_stack.append(cstack[np.where(working_true)[0]])
        
    
    def decay_particles(self):
        self.working_stack.clear()
        self.number_of_decays += self.decay_driver.run_decay(self.decay_stack, self.working_stack)
        self.decay_stack.clear()

        
        
    def get_decaying_particles(self):
        return self.decay_stack
    
    def get_final_particles(self):
        return self.final_stack
    
    
if __name__ == "__main__":
    cas_driver = CascadeDriver()
    cas_driver.run(pdg = 2212, energy = 1e5)
    
    # cas_driver.working_stack.push(pid = 2212)
    # cas_driver.working_stack.push(pid = 11, energy = 789, xdepth = 0)
    # cas_driver.working_stack.push(pid = np.array([-11, 342, 13]), energy = np.array([20, 20, 20]))
    
    # cas_driver.filter_final_particles()
    # cas_driver.filter_decaying_particles()
    
    # print(cas_driver.working_stack.valid().pid)
    # print(cas_driver.get_final_particles().valid().pid)
    # print(cas_driver.get_decaying_particles().valid().pid)
    # print(cas_driver.working_stack.valid().pid)
    # print(cas_driver.working_stack.valid().valid_code)
    
    # a = np.array([0, 0, 0], dtype=np.int64)
    # a[(0,1),] = 9
    # print(a[[0, 2]])
    
    # class Bb:
    #     def __init__(self):
    #         a = None
    
    # b = Bb()
    # setattr(b, "a", a[0:2])
    # b.a.fill(6)
    # print(a)
    

        
        
    
    
    
    