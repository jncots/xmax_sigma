import numpy as np
from particle_array import ParticleArray, FilterCode
from pdg_pid_map import PdgLists
from particle_xdepths import DefaultXdepthGetter
from hadron_inter import HadronInteraction
from decay_driver import DecayDriver
import numpy as np
import time


class CascadeDriver:
    def __init__(self, xdepth_getter = DefaultXdepthGetter(), 
                 threshold_energy = 1e3,
                 mceq_decaying_pdgs = [111]):
        
           
        self.threshold_energy = threshold_energy
        self.xdepth_getter = xdepth_getter
        self.hadron_interaction = HadronInteraction()
        self.decay_driver = DecayDriver(self.xdepth_getter)
        # Should be called after setting self.decay_driver
        self.set_decaying_pdgs(mceq_decaying_pdgs)
        
        
        self.working_stack = ParticleArray()
        self.final_stack = ParticleArray()
        self.decay_stack = ParticleArray()
        self.inter_stack = ParticleArray()
        self.failed_stack = ParticleArray()
        self.children_stack = ParticleArray()
    
    
    def set_decaying_pdgs(self, mceq_decaying_pdgs):
        pdg_lists = PdgLists()
        self.mceq_decaying_pdgs = np.array([mceq_decaying_pdgs], dtype = np.int32)
        
        # Make all mceq particles stable except the ones set as decaying
        self.stable_pdgs = (pdg_lists.mceq_particles[np.where(np.logical_not(
                            np.isin(pdg_lists.mceq_particles, self.mceq_decaying_pdgs)))[0]])
        
        # Make all pdgs decaying except the ones defined as stable
        self.decaying_pdgs = (pdg_lists.pdgs_below_abs6000[np.where(np.logical_not(
                            np.isin(pdg_lists.pdgs_below_abs6000, self.stable_pdgs)))[0]])
        
        
        self.decay_driver.set_decaying_pdgs(self.decaying_pdgs)
    
    
    def start_accumulate(self):
        self.initial_run = True
        self.accumulate_runs = True   
            
    def stop_accumulate(self):
        self.initial_run = True
        self.accumulate_runs = False
                
    
    def run(self, *, pdg, energy, 
            xdepth = 0, 
            mceq_decaying_pdgs = None, 
            threshold_energy = None):
        
        
        if threshold_energy is not None:
            self.threshold_energy = threshold_energy
        self.initial_energy = energy
        self.initial_pdg = pdg
        
        if mceq_decaying_pdgs is not None:
            self.set_decaying_pdgs(mceq_decaying_pdgs)
        
        self.number_of_decays = 0
        self.number_of_interactions = 0
        
        if self.initial_run:
            self.final_stack.clear()
            self.loop_execution_time = 0
            self.runs_number = 0
            
            if self.accumulate_runs:
                self.initial_run = False  
            
        
        self.decay_stack.clear()
        self.working_stack.clear()
        
        self.working_stack.push(pid = pdg, 
                         energy = energy, 
                         xdepth = xdepth,
                         generation_num = 0)
        
        iloop = 1
        
        
        start_time = time.time()        
        while len(self.working_stack) > 0:
            
            print(f"{iloop} Number of inter = {self.number_of_interactions}"
                  f" number of decays = {self.number_of_decays}")
            
            # self.search_for_muons(1)
            self.set_xdepths()   
            
            # self.search_for_muons(2)         
            self.run_interaction()
            
            # self.search_for_muons(3)
            self.filter_particles()
            
            # self.search_for_muons(4)
            if len(self.working_stack) == 0:
                self.decay_particles()
            
            iloop += 1
        
        self.loop_execution_time += time.time() - start_time
        self.runs_number += 1
    
    def contains_muons(self,stack, name):
        muon_number = len(np.where(np.isin(stack.valid().pid, [-13, 13]))[0])
        if muon_number > 1:
            print(f"{name} has {muon_number} muons")
        
    
    def search_for_muons(self, number):
        
        self.contains_muons(self.working_stack, f"{number} working_stack")
        self.contains_muons(self.final_stack, f"{number} final_stack")
        self.contains_muons(self.decay_stack, f"{number} decay_stack")
        self.contains_muons(self.inter_stack, f"{number} inter_stack")
        self.contains_muons(self.failed_stack, f"{number} failed_stack")
        self.contains_muons(self.children_stack, f"{number} children_stack")

        
            
    
    
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
        
        
        # Sort particles at the surface
        at_surface = np.logical_and(pvalid.xdepth_inter >= max_xdepth, pvalid.xdepth_decay >= max_xdepth)
        
        should_decay = np.isin(pvalid.pid, self.decaying_pdgs)
        should_not_decay = np.logical_not(should_decay)
        
        to_decay_stack = np.logical_and(should_decay, at_surface)
        to_final_stack = np.logical_and(should_not_decay, at_surface)
        
        self.decay_stack.append(pvalid[np.where(to_decay_stack)[0]])        
        self.final_stack.append(pvalid[np.where(to_final_stack)[0]])

        # Sort particles which are still in the atmosphere
        self.inter_stack = pvalid[np.where(pvalid.xdepth_inter < pvalid.xdepth_decay)[0]]
        self.decay_stack.append(pvalid[np.where(pvalid.xdepth_decay < pvalid.xdepth_inter)[0]])
        
        
        
    
    def run_interaction(self):
        self.children_stack.clear()
        self.failed_stack.clear()
        self.number_of_interactions += self.hadron_interaction.run_event_generator(
                                                    parents = self.inter_stack, 
                                                    children = self.children_stack, 
                                                    failed_parents = self.failed_stack)
        
              
        if len(self.failed_stack) > 0:
            fvalid = self.failed_stack.valid()
            should_decay = np.isin(fvalid.pid, self.decaying_pdgs)
            should_not_decay = np.logical_not(should_decay)
                    
            self.final_stack.append(fvalid[np.where(should_not_decay)])
            decaying_particles = fvalid[np.where(should_decay)]
            decaying_particles.filter_code[:] = FilterCode.XD_DECAY_OFF.value
            self.decay_stack.append(decaying_particles)
            
           
         
    def filter_particles(self):        
        cstack = self.children_stack.valid()
        
        below_threshold = cstack.energy < self.threshold_energy
        above_threshold = np.logical_not(below_threshold)
        
        should_decay = np.isin(cstack.pid, self.decaying_pdgs)
        should_not_decay = np.logical_not(should_decay)
        
        
        to_decay_stack = np.logical_and(below_threshold, should_decay)
        to_final_stack = np.logical_and(below_threshold, should_not_decay)
        self.decay_stack.append(cstack[np.where(to_decay_stack)[0]])
        self.final_stack.append(cstack[np.where(to_final_stack)[0]])
        
        self.working_stack.clear()
        self.working_stack.append(cstack[np.where(above_threshold)[0]])
        
        
    
    def decay_particles(self):
        if len(self.decay_stack) == 0:
            return
        self.working_stack.clear()
        self.number_of_decays += self.decay_driver.run_decay(self.decay_stack, 
                                                             final_particles=self.working_stack,
                                                             stable_particles=self.final_stack)
        self.decay_stack.clear()

        
        
    def get_decaying_particles(self):
        return self.decay_stack
    
    def get_final_particles(self):
        return self.final_stack
    
    
if __name__ == "__main__":
    cas_driver = CascadeDriver(threshold_energy = 1e3)
    cas_driver.run(pdg = 2212, energy = 1e7)

    

        
        
    
    
    
    