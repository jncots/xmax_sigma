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
        
        self.final_stack = ParticleArray()
        self.archival_stack = ParticleArray()
        
        self.working_stack = ParticleArray()
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
        
        iloop = 1
        
        
        start_time = time.time()        
        while len(self.working_stack) > 0:
            
            print(f"{iloop} Number of inter = {self.number_of_interactions}"
                  f" number of decays = {self.number_of_decays}")
            
            self.group_particles_by_slant_depth()         
            self.launch_hadron_interactions()
            self.filter_particles()
            if len(self.working_stack) == 0:
                self.launch_particle_decay()
            
            iloop += 1
        
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
        assert (init_particles == final_particles, 
                f"{init_particles - final_particles} missing particles!!!")    
                
        
            
            
    
    
    def group_particles_by_slant_depth(self):
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
        # max_xdepth = self.xdepth_getter.max_xdepth
        max_xdepth = self.stop_xdepth
        
        # stop_xdepth_arr = np.empty_like(pvalid.xdepth_inter)
        # stop_xdepth_arr.fill(self.stop_xdepth)
        
        # comparison_array = np.array([pvalid.xdepth_inter, pvalid.xdepth_decay, stop_xdepth_arr])
        # pvalid.xdepth_stop[:] = comparison_array.min(0)
        # selection_code = comparison_array.argmin(0)
        
        # self.inter_stack = pvalid[np.where(selection_code == 0)[0]]
        # self.decay_stack.append(pvalid[np.where(selection_code == 1)[0]])
        
        # hit_stop_stack = pvalid[np.where(selection_code == 2)[0]]
        # should_decay = np.isin(hit_stop_stack.pid, self.force_decaying_pdgs)
        # should_not_decay = np.logical_not(should_decay)
        
        # self.final_stack.append(hit_stop_stack[np.where(should_not_decay)[0]])
        # self.decay_stack.append(hit_stop_stack[np.where(should_decay)[0]])
        
        
        # Sort particles at the surface
        at_surface = np.logical_and(pvalid.xdepth_inter >= max_xdepth, 
                                    pvalid.xdepth_decay >= max_xdepth)
        not_at_surface = np.logical_not(at_surface)
        
        should_decay = np.isin(pvalid.pid, self.force_decaying_pdgs)
        should_not_decay = np.logical_not(should_decay)
        
        to_decay_stack_s = np.logical_and(should_decay, at_surface)
        to_final_stack = np.logical_and(should_not_decay, at_surface)
        
        
        dstack_portion_s = pvalid[np.where(to_decay_stack_s)[0]]
        dstack_portion_s.xdepth_stop[:] = max_xdepth
        self.decay_stack.append(dstack_portion_s)
        
        fstack_portion = pvalid[np.where(to_final_stack)[0]]
        fstack_portion.xdepth_stop[:] = max_xdepth
        # Hit surface, final code == 1
        fstack_portion.valid().final_code[:] = 1
        
        # print(f"Muons that hit surface = {len(np.where(fstack_portion)[0])}")        
        self.final_stack.append(fstack_portion)

        # Sort particles which are still in the atmosphere
        
        to_inter_stack = np.logical_and(not_at_surface, pvalid.xdepth_inter < pvalid.xdepth_decay)
        to_decay_stack = np.logical_and(not_at_surface, pvalid.xdepth_decay < pvalid.xdepth_inter)
        
        istack_portion = pvalid[np.where(to_inter_stack)[0]]
        istack_portion.xdepth_stop[:] = istack_portion.xdepth_inter
        self.inter_stack = istack_portion
        
        dstack_portion = pvalid[np.where(to_decay_stack)[0]]
        dstack_portion.xdepth_stop[:] = dstack_portion.xdepth_decay
        self.decay_stack.append(dstack_portion)
        
        self.search_missing_particles(init_stacks = [pvalid], 
                                      final_stacks = [fstack_portion, 
                                                      dstack_portion_s,
                                                      istack_portion,
                                                      dstack_portion,])
        
        
        
    
    def launch_hadron_interactions(self):
        if len(self.inter_stack) == 0:
            return
        
        self.children_stack.clear()
        self.rejection_stack.clear()
        self.number_of_interactions += self.hadron_interaction.run_event_generator(
                                                    parents = self.inter_stack, 
                                                    children = self.children_stack, 
                                                    failed_parents = self.rejection_stack)
        
        
        self.id_generator.generate_ids(self.children_stack.valid().id)
        self.children_stack.valid().production_code[:] = 2
        
        # Filter particles participated in interactions
        parents_true = np.logical_not(np.isin(self.inter_stack.valid().pid, 
                                              self.rejection_stack.valid().pid))
        parents = self.inter_stack[np.where(parents_true)[0]]
        parents.valid().final_code[:] = 2
        # And record them in archival stack
        self.archival_stack.append(parents)
        

        parent_energy = np.sum(self.inter_stack.valid().energy)
        failed_energy = np.sum(self.rejection_stack.valid().energy)
        children_energy = np.sum(self.children_stack.valid().energy)
        secondary_energy = children_energy + failed_energy
        print(f"Parent energy = {parent_energy}")
        print(f"Failed energy = {failed_energy}")  
        print(f"Children = {children_energy}")
        print(f"Secondary energy = {secondary_energy}") 
        print(f"Energy conservation = {100*(secondary_energy - parent_energy)/parent_energy} %")
              
        if len(self.rejection_stack) > 0:
            fvalid = self.rejection_stack.valid()
            should_decay = np.isin(fvalid.pid, self.force_decaying_pdgs)
            should_not_decay = np.logical_not(should_decay)
            
            
            fstack_portion = fvalid[np.where(should_not_decay)]
            # fstack_portion.xdepth_stop[:] = self.stop_xdepth        
            self.final_stack.append(fstack_portion)
            
            decaying_particles = fvalid[np.where(should_decay)]
            # decaying_particles.filter_code[:] = FilterCode.XD_DECAY_OFF.value
            self.decay_stack.append(decaying_particles)
        
            self.search_missing_particles(init_stacks = [self.rejection_stack.valid()], 
                                        final_stacks = [fstack_portion, 
                                                        decaying_particles,
                                                        ])   
           
         
    def filter_particles(self):        
        cstack = self.children_stack.valid()
        
        below_threshold = cstack.energy < self.threshold_energy
        above_threshold = np.logical_not(below_threshold)
        
        should_decay = np.isin(cstack.pid, self.force_decaying_pdgs)
        should_not_decay = np.logical_not(should_decay)
        
        
        to_decay_stack = np.logical_and(below_threshold, should_decay)
        to_final_stack = np.logical_and(below_threshold, should_not_decay)
        
        decay_portion = cstack[np.where(to_decay_stack)[0]]
        self.decay_stack.append(decay_portion)
        
        fstack_portion = cstack[np.where(to_final_stack)[0]]
        fstack_portion.xdepth_stop[:] = self.stop_xdepth
        self.final_stack.append(fstack_portion)
        
        self.working_stack.clear()
        
        work_portion = cstack[np.where(above_threshold)[0]] 
        self.working_stack.append(work_portion)
        
        self.search_missing_particles(init_stacks = [cstack], 
                                        final_stacks = [fstack_portion,
                                                        decay_portion, 
                                                        work_portion,
                                                        ])   
        
        
    
    def launch_particle_decay(self):
        if len(self.decay_stack) == 0:
            return
        
        self.rejection_stack.clear()
        self.working_stack.clear()
        self.number_of_decays += self.decay_driver.run_decay(self.decay_stack, 
                                                             decayed_particles=self.working_stack,
                                                             stable_particles=self.rejection_stack)
        
        
        print(f"Decaying stack = {self.decay_stack.valid().pid}")
        print(f"Decayed particles = {self.working_stack.valid().pid}")
        print(f"Stable parts = {self.rejection_stack.valid().pid}")
        
        init_energy = np.sum(self.decay_stack.valid().energy)
        fin_dec_energy = np.sum(self.working_stack.valid().energy)
        stab_energy = np.sum(self.rejection_stack.valid().energy)
        
        print(f"Init = {init_energy}, fin_dec = {fin_dec_energy}, stab = {stab_energy}")
        print(f"Conservation of energy = {((fin_dec_energy + stab_energy) -  init_energy)/init_energy} %")
        
        
        # Filter particles participated in decay
        parents_true = np.logical_not(np.isin(self.decay_stack.valid().pid, 
                                              self.rejection_stack.valid().pid))
        parents = self.decay_stack[np.where(parents_true)[0]]
        # Final_code = 3 means decay
        parents.valid().final_code[:] = 3
        # And record them in archival stack
        self.archival_stack.append(parents)
        
        self.id_generator.generate_ids(self.working_stack.valid().id)
        self.rejection_stack.valid().xdepth_stop[:] = self.stop_xdepth
        self.final_stack.append(self.rejection_stack)
        
        
        # self.search_missing_particles(init_stacks = [self.rejection_stack.valid()], 
        #                                 final_stacks = [fstack_portion,
        #                                                 decay_portion, 
        #                                                 work_portion,
        #                                                 ]) 
        
        
        self.decay_stack.clear()

        
        
    def get_decaying_particles(self):
        return self.decay_stack
    
    def get_final_particles(self):
        return self.final_stack        
    
if __name__ == "__main__":
    cas_driver = CascadeDriver(threshold_energy = 1e3)
    cas_driver.run(pdg = 2212, energy = 1e7)

    

        
        
    
    
    
    