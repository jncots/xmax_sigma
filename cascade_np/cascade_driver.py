import numpy as np
from particle_array import ParticleArray, FilterCode
from pdg_pid_map import PdgLists
from particle_xdepths import DefaultXdepthGetter
from hadron_inter import HadronInteraction
from decay_driver import DecayDriver
import numpy as np
import time
from id_generator import IdGenerator


class CascadeDriver:
    def __init__(self, zenith_angle):
        
        self.id_generator = IdGenerator()
        self.xdepth_getter = DefaultXdepthGetter(zenith_angle)
        self.hadron_interaction = HadronInteraction()
        self.decay_driver = DecayDriver(self.xdepth_getter)
        
        
        self.final_stack_decay = ParticleArray()
        self.final_stack = ParticleArray()
        self.archival_stack = ParticleArray()
        self.generated_stack = ParticleArray()
        
        self.working_stack = ParticleArray()
        self.working_stack_filter = ParticleArray()
        self.decay_stack = ParticleArray()
        self.inter_stack = ParticleArray()
        self.rejection_stack = ParticleArray()
        self.children_stack = ParticleArray()
    
    
        
    def set_decaying_pdgs(self, mceq_decaying_pdgs):
        pdg_lists = PdgLists()
        
        # Force decaying pdgs are particles which do not present in mceq
        # and should decay event if they hit final conditions
        # Every particle which is not in mceq
        self.force_decaying_pdgs = (pdg_lists.pdgs_below_abs6000[np.where(np.logical_not(
                            np.isin(pdg_lists.pdgs_below_abs6000, pdg_lists.mceq_particles)))[0]])
        
        # Natural decaying pdgs are particles which present in mceq
        # and can be present in final stack
        # Mceq particles that should decay
        self.natural_decaying_pdgs = np.array([mceq_decaying_pdgs], dtype = np.int32)
        
        
        # Mceq particles that is stable
        self.stable_pdgs = (pdg_lists.mceq_particles[np.where(np.logical_not(
                            np.isin(pdg_lists.mceq_particles, self.natural_decaying_pdgs)))[0]])
        
        # Make all pdgs decaying except the ones defined as stable
        decaying_pdgs = (pdg_lists.pdgs_below_abs6000[np.where(np.logical_not(
                            np.isin(pdg_lists.pdgs_below_abs6000, self.stable_pdgs)))[0]])
        
        
        self.decay_driver.set_decaying_pdgs(decaying_pdgs=decaying_pdgs,
                                            stable_pdgs=self.stable_pdgs)
    
    
    
    def simulation_parameters(self, *, pdg, energy, 
                              threshold_energy,
                              mceq_decaying_pdgs,
                              xdepth = 0,
                              stop_height = 0,
                              accumulate_runs = False,
                              reset_ids = False):
            
        self.initial_pdg = pdg
        self.initial_energy = energy
        self.threshold_energy = threshold_energy
        self.initial_xdepth = xdepth
        
        self.set_decaying_pdgs(mceq_decaying_pdgs)
        
        self.stop_xdepth = self.xdepth_getter.xdepth_conversion.convert_h2x(stop_height * 1e5)
        self.xdepth_getter.set_stop_xdepth(self.stop_xdepth)
        print(f"stop depth = {self.stop_xdepth}")
        
        if reset_ids:
            self.id_generator = IdGenerator()
        
        self.initial_run = True    
        self.accumulate_runs = accumulate_runs           
    
    def run(self):
        
        self.working_stack.clear()
        self.decay_stack.clear()
                
        if self.initial_run:
            self.final_stack.clear()
            self.archival_stack.clear()
            self.generated_stack.clear()
            self.number_of_decays = 0
            self.number_of_interactions = 0
            self.loop_execution_time = 0
            self.runs_number = 0
            
            if self.accumulate_runs:
                self.initial_run = False  
            
        
        self.working_stack.push(pid = self.initial_pdg, 
                         energy = self.initial_energy, 
                         xdepth = self.initial_xdepth,
                         generation_num = 0)
        
        self.id_generator.generate_ids(self.working_stack.valid().id)
        self.generated_stack.append(self.working_stack)
        
        iloop = 1
        
        
        start_time = time.time()        
        while len(self.working_stack) > 0:
            
            print(f"\r{iloop} Number of inter = {self.number_of_interactions}"
                  f" number of decays = {self.number_of_decays}")
            
            self.filter_by_energy()
            self.filter_by_slant_depth()         
            self.run_hadron_interactions()
            if len(self.working_stack) == 0:
                self.run_particle_decay()
            
            iloop += 1
        
        # self.run_decay_at_surface()
        
        self.loop_execution_time += time.time() - start_time
        self.runs_number += 1
    
    def contains_muons(self, stack, name):
        muon_number = len(np.where(np.isin(stack.valid().pid, [-13, 13]))[0])
        if muon_number > 1:
            print(f"{name} has {muon_number} muons")
        
    
    def search_for_muons(self, number):
        
        self.contains_muons(self.working_stack, f"{number} working_stack")
        self.contains_muons(self.final_stack, f"{number} final_stack")
        self.contains_muons(self.decay_stack, f"{number} decay_stack")
        self.contains_muons(self.inter_stack, f"{number} inter_stack")
        self.contains_muons(self.rejection_stack, f"{number} rejection_stack")
        self.contains_muons(self.children_stack, f"{number} children_stack")
        
    def final_contains_small_xdepth_stop(self, stack, name):
        entry_number = len(np.where(stack.valid().xdepth_stop < self.stop_xdepth)[0])
        if entry_number > 0:
            print(f"{name} has {entry_number} entries")
        
    
    def search_for_xdepth_stop(self, number):
        
        self.final_contains_small_xdepth_stop(self.final_stack, f"{number} final_stack")
        # self.contains_muons(self.decay_stack, f"{number} decay_stack")
        # self.contains_muons(self.inter_stack, f"{number} inter_stack")
        # self.contains_muons(self.rejection_stack, f"{number} rejection_stack")
        # self.contains_muons(self.children_stack, f"{number} children_stack")    

    
    def search_missing_particles(self, init_stacks, final_stacks):
        
        init_particles = 0
        for st in init_stacks:
            init_particles += len(st)
        
        final_particles = 0
        for st in final_stacks:
            final_particles += len(st)
        
        
        print(f"TEST: {init_particles - final_particles} MISSING,"
              f" INIT = {init_particles}, FINAL = {final_particles}")    
        assert init_particles == final_particles, f"{init_particles - final_particles} missing particles!!!"    
                
        
    
    def filter_by_energy(self):   
        
        wstack = self.working_stack.valid()
        
        above_threshold = wstack.energy > self.threshold_energy
        below_threshold =  np.logical_not(above_threshold)
        
        above_threshold = np.where(above_threshold)[0]
        self.working_stack_filter.clear()
        self.working_stack_filter.append(wstack[above_threshold])
        
        
        should_decay = np.isin(wstack.pid, self.force_decaying_pdgs)
        should_not_decay = np.logical_not(should_decay)
        
        should_decay = np.where(np.logical_and(below_threshold, should_decay))[0]
        should_not_decay = np.where(np.logical_and(below_threshold, should_not_decay))[0]

        self.final_stack_decay.append(wstack[should_decay])
        self.final_stack.append(wstack[should_not_decay])     
        
        
        particle_number = above_threshold.size + should_decay.size + should_not_decay.size
        assert particle_number == len(wstack), (
                "Number of distributed particles not equal to number of initial particles")
        
        self.working_stack.clear()
    
    
    def filter_by_slant_depth(self):
        """pstack is assumed to contain particles which:
            * has energy > threshold_energy
            * xdepth < surface_xdepth
            * has pdg which can be can decay or interact

        Args:
            pstack (_type_): _description_

        Returns:
            _type_: _description_
        """
        
        wstack = self.working_stack_filter.valid()
        self.xdepth_getter.get_decay_xdepth(wstack)
        self.xdepth_getter.get_inter_xdepth(wstack)
        
        max_xdepth = self.stop_xdepth
                
        # Sort particles at the surface
        at_surface = np.logical_and(wstack.xdepth_inter >= max_xdepth, 
                                    wstack.xdepth_decay >= max_xdepth)
        
        should_decay = np.isin(wstack.pid, self.force_decaying_pdgs)
        should_not_decay = np.logical_not(should_decay)
        
        should_decay = np.where(np.logical_and(at_surface, should_decay))[0]
        should_not_decay = np.where(np.logical_and(at_surface, should_not_decay))[0]
        
        
        # Particles that should be decayed at surface
        dstack_portion = wstack[should_decay]
        dstack_portion.xdepth_stop[:] = max_xdepth
        dstack_portion.final_code[:] = 1
        self.final_stack_decay.append(dstack_portion)
        
        # Particles that are already at their final stage
        fstack_portion = wstack[should_not_decay]
        fstack_portion.xdepth_stop[:] = max_xdepth
        fstack_portion.final_code[:] = 1
        self.final_stack.append(fstack_portion)

        # Sort particles which are still in the atmosphere
        istack_portion = wstack[wstack.xdepth_inter < wstack.xdepth_decay]
        istack_portion.xdepth_stop[:] = istack_portion.xdepth_inter
        self.inter_stack = istack_portion
        
        dstack_portion = wstack[wstack.xdepth_decay > wstack.xdepth_inter]
        dstack_portion.xdepth_stop[:] = dstack_portion.xdepth_decay
        dstack_portion.final_code[:] = -1
        self.decay_stack.append(dstack_portion)
        
        
        
    
    def run_hadron_interactions(self):
        self.children_stack.clear()
        self.rejection_stack.clear()
        
        if len(self.inter_stack) == 0:
            return
        
        self.number_of_interactions += self.hadron_interaction.run_event_generator(
                                                    parents = self.inter_stack, 
                                                    children = self.children_stack, 
                                                    failed_parents = self.rejection_stack)
        
        
        self.id_generator.generate_ids(self.children_stack.valid().id)
        self.children_stack.valid().production_code[:] = 2
        
        self.generated_stack.append(self.children_stack)
        
        # Filter particles participated in interactions
        parents_true = np.logical_not(np.isin(self.inter_stack.valid().id, 
                                              self.rejection_stack.valid().id))
        parents = self.inter_stack[np.where(parents_true)[0]]
        parents.valid().final_code[:] = 2
        # And record them in archival stack
        self.archival_stack.append(parents)
        
        if len(self.rejection_stack) > 0:
            self.archival_stack.append(self.rejection_stack)
        

        self.working_stack.clear()
        self.working_stack.append(self.children_stack)
              
        # if len(self.rejection_stack) > 0:
        #     fvalid = self.rejection_stack.valid()
        #     should_decay = np.isin(fvalid.pid, self.force_decaying_pdgs)
        #     should_not_decay = np.logical_not(should_decay)
            
            
        #     fstack_portion = fvalid[np.where(should_not_decay)]
        #     # fstack_portion.xdepth_stop[:] = self.stop_xdepth  
        #     # self.debug_append_final_stack(fstack_portion)       
        #     self.final_stack.append(fstack_portion)
            
        #     decaying_particles = fvalid[np.where(should_decay)]
        #     # decaying_particles.filter_code[:] = FilterCode.XD_DECAY_OFF.value
        #     self.decay_stack.append(decaying_particles)
        
        #     # self.search_missing_particles(init_stacks = [self.rejection_stack.valid()], 
        #     #                             final_stacks = [fstack_portion, 
        #     #                                             decaying_particles,
        #     #                                             ])     
         
    # def filter_particles(self):   
    #     cstack = self.children_stack.valid()
        
    #     below_threshold = cstack.energy < self.threshold_energy
    #     above_threshold = np.logical_not(below_threshold)
        
    #     should_decay = np.isin(cstack.pid, self.force_decaying_pdgs)
    #     should_not_decay = np.logical_not(should_decay)
        
        
    #     to_decay_stack = np.logical_and(below_threshold, should_decay)
    #     to_final_stack = np.logical_and(below_threshold, should_not_decay)
        
    #     decay_portion = cstack[np.where(to_decay_stack)[0]]
    #     decay_portion.final_code[:] = 4
    #     self.decay_stack.append(decay_portion)
        
    #     fstack_portion = cstack[np.where(to_final_stack)[0]]
    #     fstack_portion.xdepth_stop[:] = self.stop_xdepth
    #     # self.debug_append_final_stack(fstack_portion) 
    #     self.final_stack.append(fstack_portion)
        
    #     self.working_stack.clear()
        
    #     work_portion = cstack[np.where(above_threshold)[0]] 
    #     self.working_stack.append(work_portion)
        
        # print(f"cstack.id = {cstack.valid().id} ,cstack.pid = {cstack.valid().pid}, cstack.energy = {cstack.valid().energy}")
        # print(f"fstack_portion.pid = {fstack_portion.valid().pid}, fstack_portion.energy = {fstack_portion.valid().energy}")
        # print(f"decay_portion.pid = {decay_portion.valid().pid}, decay_portion.energy = {decay_portion.valid().energy}")
        # print(f"work_portion.id = {work_portion.valid().id} ,work_portion.pid = {work_portion.valid().pid}, work_portion.energy = {work_portion.valid().energy}")
        
        # self.search_missing_particles(init_stacks = [cstack], 
        #                                 final_stacks = [fstack_portion,
        #                                                 decay_portion, 
        #                                                 work_portion,
        #                                                 ])   
    
    def run_decay_at_surface(self):
        if len(self.final_stack_decay) == 0:
            return
        
        self.rejection_stack.clear()
        self.working_stack.clear()
        
        self.number_of_decays += self.decay_driver.run_decay(self.final_stack_decay, 
                                                             decayed_particles=self.working_stack,
                                                             stable_particles=self.rejection_stack)
        
        
        self.final_stack.append(self.working_stack)
        self.final_stack.append(self.rejection_stack)
        
        
        
    
    def run_particle_decay(self):
        if len(self.decay_stack) == 0:
            return
        
        self.rejection_stack.clear()
        self.working_stack.clear()
        self.number_of_decays += self.decay_driver.run_decay(self.decay_stack, 
                                                             decayed_particles=self.working_stack,
                                                             stable_particles=self.rejection_stack)
        
        
        # Filter particles participated in decay
        parents_true = np.logical_not(np.isin(self.decay_stack.valid().id, 
                                              self.rejection_stack.valid().id))
        parents = self.decay_stack[np.where(parents_true)[0]]
        # Final_code = 3 means decay
        parents.valid().final_code[:] = 3
        # And record them in archival stack
        self.archival_stack.append(parents)
        
        self.id_generator.generate_ids(self.working_stack.valid().id)
        self.generated_stack.append(self.working_stack)
        
        self.rejection_stack.valid().xdepth_stop[:] = self.stop_xdepth
        
        
        # self.debug_append_final_stack(self.rejection_stack) 
        self.final_stack.append(self.rejection_stack) 
        
        # Use rejection_stack again
        self.rejection_stack.clear()
        should_be_finals = np.isin(self.working_stack.valid().final_code, np.array([1, 4], dtype=np.int32))
        should_be_processed = np.logical_not(should_be_finals)
        self.final_stack.append(self.working_stack[np.where(should_be_finals)[0]])
        self.rejection_stack.append(self.working_stack[np.where(should_be_processed)[0]])
        self.working_stack.clear()
        
        self.working_stack.append(self.rejection_stack)  
        
        self.decay_stack.clear()

    def debug_append_final_stack(self, array_to_append):
        print(f"START DEBUG_APPEND")
        print(f"final.id = {self.final_stack.valid().id}")
        
        cont = np.where(np.isin(array_to_append.valid().id, self.final_stack.valid().id))[0]
        if len(cont) > 0:
            print(f"IDs = {array_to_append.valid().id[cont]},\n PDGS = {array_to_append.valid().pid[cont]}"
                  f" ENs = {array_to_append.valid().energy[cont]}")
            print(f"final.id = {self.final_stack.valid().id}")
            raise ValueError("WRONG APPEND")
        
        print(f"END DEBUG_APPEND")
        
    def get_decaying_particles(self):
        return self.decay_stack
    
    def get_final_particles(self):
        return self.final_stack        
    
if __name__ == "__main__":
    cas_driver = CascadeDriver(threshold_energy = 1e3)
    cas_driver.run(pdg = 2212, energy = 1e7)

    

        
        
    
    
    
    