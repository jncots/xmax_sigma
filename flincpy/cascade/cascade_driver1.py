import numpy as np

from data_structs.particle_array import ParticleArray
from data_structs.pdg_pid_map import PdgLists
from data_structs.id_generator import IdGenerator

from propagation.particle_xdepths import DefaultXdepthGetter

from process.hadron_inter import HadronInteraction
from process.decay_driver import DecayDriver
import numpy as np
import time



# def example_interaction_model():
#     import chromo
#     target = chromo.kinematics.CompositeTarget([("N", 0.78), ("O", 0.22)])
#     ekin = chromo.kinematics.FixedTarget(1e7, "O16", target)
#     model = chromo.models.DpmjetIII191
    
#     return InteractionModel(model, ekin, target)

class InteractionModel:
    def __init__(self, model, initial_kinematics, target):
        self.target = target
        self.event_generator = model(initial_kinematics)

class CascadeDriver:
    def __init__(self, imodel):
        self.imodel = imodel
        self.id_generator = IdGenerator()
        self.pdg_lists = PdgLists()
        
        
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
    
    
    def set_zenith_angle(self, zenith_angle):
        self.xdepth_getter = DefaultXdepthGetter(zenith_angle)
        self.decay_driver = DecayDriver(self.xdepth_getter)
        self.hadron_interaction = HadronInteraction(self.imodel, self.xdepth_getter)
        
        
    def set_decaying_pdgs(self, mceq_decaying_pdgs):
        
        # Force decaying pdgs are particles which do not present in mceq
        # and should decay event if they hit final conditions
        # Every particle which is not in mceq
        self.force_decaying_pdgs = (self.pdg_lists.pdgs_below_abs6000[np.where(np.logical_not(
                            np.isin(self.pdg_lists.pdgs_below_abs6000, self.pdg_lists.mceq_particles)))[0]])
        
        # Natural decaying pdgs are particles which present in mceq
        # and can be present in final stack
        # Mceq particles that should decay
        self.natural_decaying_pdgs = np.array([mceq_decaying_pdgs], dtype = np.int32)
        
        
        # Mceq particles that is stable
        self.stable_pdgs = (self.pdg_lists.mceq_particles[np.where(np.logical_not(
                            np.isin(self.pdg_lists.mceq_particles, self.natural_decaying_pdgs)))[0]])
        
        # Make all pdgs decaying except the ones defined as stable
        decaying_pdgs = (self.pdg_lists.pdgs_below_abs6000[np.where(np.logical_not(
                            np.isin(self.pdg_lists.pdgs_below_abs6000, self.stable_pdgs)))[0]])
        
        
        self.decay_driver.set_decaying_pdgs(decaying_pdgs=decaying_pdgs,
                                            stable_pdgs=self.stable_pdgs)
    
    
    
    def simulation_parameters(self, *, pdg, energy, 
                              threshold_energy,
                              zenith_angle,
                              mceq_decaying_pdgs,
                              xdepth = 0,
                              stop_height = 0,
                              accumulate_runs = False,
                              reset_ids = False):
            
        self.initial_pdg = pdg
        self.initial_energy = energy
        self.threshold_energy = threshold_energy
        self.initial_xdepth = xdepth
        
        self.set_zenith_angle(zenith_angle)
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
            self.final_stack_decay.clear()
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
        
        self.run_decay_at_surface()
        
        self.loop_execution_time += time.time() - start_time
        self.runs_number += 1
    
    
    def filter_by_energy(self):   
        
        wstack = self.working_stack.valid()
        
        # Put in the working stack particles with above threshold energy
        is_above_threshold = wstack.energy > self.threshold_energy
        above_threshold_idx = np.where(is_above_threshold)[0]
        self.working_stack_filter.clear()
        self.working_stack_filter.append(wstack[above_threshold_idx])
        
        
        # Separate the particles that should decay before moving to final_stack
        # from particles that are already in final state and should move to final_stack
        is_below_threshold =  np.logical_not(is_above_threshold)
        is_decaying = np.isin(wstack.pid, self.force_decaying_pdgs)
        is_stable = np.logical_not(is_decaying)
        
        is_decaying_idx = np.where(np.logical_and(is_below_threshold, is_decaying))[0]
        is_stable_idx = np.where(np.logical_and(is_below_threshold, is_stable))[0]

        self.final_stack_decay.append(wstack[is_decaying_idx])
        self.final_stack.append(wstack[is_stable_idx])     
        
        
        particle_number = above_threshold_idx.size + is_decaying_idx.size + is_stable_idx.size
        
        if above_threshold_idx.size < len(wstack):
            print(f"Above threshold = {above_threshold_idx.size/len(wstack)*100} %")
        
        assert particle_number == len(wstack), (
                "Number of distributed particles not equal to number of initial particles")
        
        self.working_stack.clear()
    
    
    def filter_by_slant_depth(self):
        wstack = self.working_stack_filter.valid()
        
        is_stable_lepton = np.isin(wstack.pid, self.pdg_lists.leptons_stable_mceq)
        is_decaying_lepton = np.isin(wstack.pid, self.pdg_lists.leptons_decay_mceq)
        is_mixing_hadron = np.isin(wstack.pid, self.pdg_lists.hadrons_mix_mceq)
        is_stable_hadron = np.isin(wstack.pid, self.pdg_lists.hadrons_stable_mceq)
        
        
        self.final_stack.append(wstack[is_stable_lepton])
        self.decay_stack.append(wstack[is_decaying_lepton])
        
        mixing_hadrons = wstack[is_mixing_hadron]
        is_mixing = self.pdg_lists.hadron_emix[mixing_hadrons.pid] < mixing_hadrons.energy
        is_only_decaying = np.logical_not(is_mixing)
        
        self.decay_stack.append(mixing_hadrons[is_only_decaying])
        
        
        self.working_stack.append(mixing_hadrons[is_mixing])
        self.working_stack.append(wstack[is_stable_hadron])
        
        wstack = self.working_stack.valid()
        
        print(f"Wstack for interaction = {wstack.pid}")
        
        
        self.xdepth_getter.get_decay_xdepth(wstack)
        self.xdepth_getter.get_inter_xdepth(wstack)
        
        max_xdepth = self.stop_xdepth
                
        # Sort particles at the surface
        at_surface = np.logical_and(wstack.xdepth_inter >= max_xdepth, 
                                    wstack.xdepth_decay >= max_xdepth)
        
        not_at_surface = np.logical_not(at_surface)
        
        should_decay = np.isin(wstack.pid, self.force_decaying_pdgs)
        should_not_decay = np.logical_not(should_decay)
        
        should_decay = np.where(np.logical_and(at_surface, should_decay))[0]
        should_not_decay = np.where(np.logical_and(at_surface, should_not_decay))[0]
        
        
        # Particles that should be decayed at surface
        dfstack_portion = wstack[should_decay]
        dfstack_portion.xdepth_stop[:] = max_xdepth
        dfstack_portion.final_code[:] = 1
        self.final_stack_decay.append(dfstack_portion)
        
        # Particles that are already at their final stage
        fstack_portion = wstack[should_not_decay]
        fstack_portion.xdepth_stop[:] = max_xdepth
        fstack_portion.final_code[:] = 1
        self.final_stack.append(fstack_portion)

        # Sort particles which are still in the atmosphere
        istack_true = np.logical_and(wstack.xdepth_inter < wstack.xdepth_decay, not_at_surface)
        istack_true = np.where(istack_true)[0]
        istack_portion = wstack[istack_true]        
        istack_portion.xdepth_stop[:] = istack_portion.xdepth_inter
        self.inter_stack = istack_portion
        
        dstack_true = np.logical_and(wstack.xdepth_inter > wstack.xdepth_decay, not_at_surface)
        dstack_true = np.where(dstack_true)[0]
        dstack_portion = wstack[dstack_true]
        dstack_portion.xdepth_stop[:] = dstack_portion.xdepth_decay
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
        
        # Take only MCEq particles
        is_mceq_particles = np.isin(self.children_stack.pid, self.pdg_lists.mceq_particles)
        self.children_stack = self.children_stack[is_mceq_particles]
        
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
            self.decay_stack.append(self.rejection_stack)
            
            rej_muons = np.where(np.isin(self.rejection_stack.valid().pid, 
                                         np.array([-13, 13], dtype = np.int32)))[0]
            if len(rej_muons) > 0:
                print(f"Muons pdgs = {self.rejection_stack.pid[rej_muons]}"
                      f" energy = {self.rejection_stack.energy[rej_muons]}"
                      f" xdepth = {self.rejection_stack.xdepth[rej_muons]}"
                      f" xdepth_decay = {self.rejection_stack.xdepth_decay[rej_muons]}"
                      f" xdepth_inter = {self.rejection_stack.xdepth_inter[rej_muons]}")
        

        self.working_stack.clear()
        self.working_stack.append(self.children_stack)

    
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
        
        
    def run_decay_at_surface(self):
        
        # print(f"Run decay1, number = {len(self.final_stack_decay)}")
        if len(self.final_stack_decay) == 0:
            return
        
        # print("Run decay2")
        self.rejection_stack.clear()
        self.working_stack.clear()
        
        self.number_of_decays += self.decay_driver.run_decay(self.final_stack_decay, 
                                                             decayed_particles=self.working_stack,
                                                             stable_particles=self.rejection_stack)
        
        
        self.final_stack.append(self.working_stack)
        self.final_stack.append(self.rejection_stack)
        
    def get_decaying_particles(self):
        return self.decay_stack
    
    def get_final_particles(self):
        return self.final_stack        
    
if __name__ == "__main__":
    cas_driver = CascadeDriver(threshold_energy = 1e3)
    cas_driver.run(pdg = 2212, energy = 1e7)

    

        
        
    
    
    
    