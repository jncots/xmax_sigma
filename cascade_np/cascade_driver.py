import numpy as np
from particle_array import ParticleArray, FilterCode
from pdg_pid_map import PdgLists
from particle_xdepths import default_xdepth_getter
from hadron_inter import HadronInteraction
import numpy as np

class CascadeDriver:
    def __init__(self, xdepth_getter = default_xdepth_getter):

        
        self.xdepth_getter = xdepth_getter
        self.pdg_lists = PdgLists()
        self.hadron_interaction = HadronInteraction()
        self.working_stack = ParticleArray(1000000)
        self.final_stack = ParticleArray(1000000)
        self.decay_stack = ParticleArray(1000000)
        self.inter_stack = ParticleArray(1000000)
    
    
    def run(self, *, pdg, energy, xdepth = 0):
        self.working_stack.push(pid = pdg, 
                         energy = energy, 
                         xdepth = xdepth)
        
        
        while len(self.working_stack) > 0:
            
            self.set_xdepths(self.working_stack.valid())
            self.run_interaction()
        #     self.cascade_event.get_event_particles(self.working_stack)
        
        # self.final_stack.append(self.working_stack)
        
    # def filter_below_threshold(self):
    #     threshold_energy = 1e3
        
    #     wstack = self.working_stack.valid()
    #     finals_slice = np.where(wstack.energy < threshold_energy & wstack.valid_code == 1)[0]
        
    #     self.decay_stack.append(wstack[finals_slice])
    #     wstack.valid_code[finals_slice] = 0
        
    
    
    def set_xdepths(self, pstack):
        """pstack is assumed to contain particles which:
            * has energy > threshold_energy
            * xdepth < surface_xdepth
            * has pdg which can be can decay or interact

        Args:
            pstack (_type_): _description_

        Returns:
            _type_: _description_
        """
        
        self.xdepth_getter.get_decay_xdepth(pstack)
        self.xdepth_getter.get_inter_xdepth(pstack)
        
        
        pst = pstack.valid()
        max_xdepth = self.xdepth_getter.max_xdepth
        
        # Copies of particles with specific properties
        self.final_stack.append(pst[np.where((pst.xdepth_inter >= max_xdepth) &
                                (pst.xdepth_decay >= max_xdepth))[0]])
        
        self.inter_stack = pst[np.where(pst.xdepth_inter < pst.xdepth_decay)[0]]
        
        

    
        
        decaying_particles = pst[np.where(pst.xdepth_inter > pst.xdepth_decay)[0]]
        # Mark that xdepth_decay is already set
        decaying_particles.filter_code[:] = FilterCode.XD_DECAY_ON.value
        self.decay_stack.append(decaying_particles)
        
    
    def run_interaction(self):
        self.hadron_interaction.run_event_generator(self.inter_stack)
           
        self.working_stack = self.hadron_interaction.get_children()
        failed = self.hadron_interaction.get_failed()
           
        if len(failed) > 0:
            in_finals = np.in1d(failed.pid, self.pdg_lists.final_pdgs)
            not_in_finals = np.logical_not(in_finals)           
            self.final_stack.append(failed[np.where(in_finals)])
            self.decay_stack.append(failed[np.where(not_in_finals)])
           
           
    
    def filter_decaying_particles(self):
        decay_pdgs = np.array([111], dtype = np.int32)
        
        wstack = self.working_stack.valid()
        finals_slice = np.where(np.in1d(wstack.pid, decay_pdgs) & (wstack.valid_code == 1))[0]
        
        self.decay_stack.append(wstack[finals_slice])
        wstack.valid_code[finals_slice] = 0    
    
    
    def filter_final_particles(self):
        threshold_energy = 1e3
        final_pdgs = np.array([-11, 11, 12, -13, 13, 14, 16, 22], dtype = np.int32)
        
        wstack = self.working_stack.valid()
        finals_slice = np.where((np.in1d(wstack.pid, self.pdg_lists.final_pdgs) 
                                | (wstack.energy < threshold_energy)) 
                                & (wstack.valid_code == 1))[0]
        
        self.final_stack.append(wstack[finals_slice])
        wstack.valid_code[finals_slice] = 0
        print("wstack", wstack.valid_code)
        
        
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
    

        
        
    
    
    
    