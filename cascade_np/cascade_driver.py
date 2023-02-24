import numpy as np
from particle_array import ParticleArray, FilterCode
from pdg_pid_map import PdgLists
from particle_xdepths import default_xdepth_getter
from hadron_inter import HadronInteraction
from decay_driver import DecayDriver
import numpy as np


def check_xdepth(stack):
    svalid = stack.valid()
    wrong = np.where(np.logical_or(svalid.xdepth < 0, svalid.xdepth > 1200))[0]
    if len(wrong):
        raise ValueError(f"wrong = {wrong}"
                         f"xdepth = {svalid.xdepth[wrong]}")

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
        self.final_stack.clear()
        self.decay_stack.clear()
        self.working_stack.clear()
        
        self.working_stack.push(pid = pdg, 
                         energy = energy, 
                         xdepth = xdepth,
                         generation_num = 0)
        
        iloop = 1
        while len(self.working_stack) > 0:
            
            if len(self.working_stack) in [1, 2, 3]:
                print("pid = ", self.working_stack.valid().pid)
                print("energy = ", self.working_stack.valid().energy)
                print("xdepth = ", self.working_stack.valid().xdepth)
                print("xdepth_decay = ", self.working_stack.valid().xdepth_decay)
                print("xdepth_inter = ", self.working_stack.valid().xdepth_inter)
                print("gen_num = ", self.working_stack.valid().generation_num)
            
            print(f"{iloop} Number of inter = {self.number_of_interactions}"
                  f" number of decays = {self.number_of_decays}")
            
            print(f"Final = {len(self.final_stack)}"
                  f" Decaying = {len(self.decay_stack)}"
                  f" Children = {len(self.children_stack)}"
                  f" Working = {len(self.working_stack)}"
                  )
            self.set_xdepths()
            check_xdepth(self.final_stack)
            check_xdepth(self.decay_stack)
            check_xdepth(self.working_stack)
            
            self.run_interaction()
            
            check_xdepth(self.final_stack)
            check_xdepth(self.decay_stack)
            check_xdepth(self.working_stack)
            
            self.filter_particles()
            
            check_xdepth(self.final_stack)
            check_xdepth(self.decay_stack)
            check_xdepth(self.working_stack)
            
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
        
        # print("xdepth = ", np.where(np.abs(self.working_stack.valid().xdepth) > 2000)[0])
        # print("pid = ", self.working_stack.pid[np.where(np.abs(self.working_stack.valid().xdepth) > 2000)[0]])
        # print("xdepth = ", self.working_stack.xdepth[np.where(np.abs(self.working_stack.valid().xdepth) > 2000)[0]])
        
        self.xdepth_getter.get_decay_xdepth(self.working_stack)
        self.xdepth_getter.get_inter_xdepth(self.working_stack)
        
        
        pvalid = self.working_stack.valid()
        max_xdepth = self.xdepth_getter.max_xdepth
        
        # Copies of particles with specific properties
        
        hit_earth = np.logical_and(pvalid.xdepth_inter >= max_xdepth,
                                   pvalid.xdepth_decay >= max_xdepth)
        
        will_interact = pvalid.xdepth_inter < pvalid.xdepth_decay
        will_decay = pvalid.xdepth_decay < pvalid.xdepth_inter
        
        res = np.where(hit_earth, 1, 0) + np.where(will_interact, 1, 0) + np.where(will_decay, 1, 0)
        
        if len(np.where(res > 1)[0]) > 0:
            raise AssertionError(f"res = {np.where(res > 1)[0]}\n"
                                 f" pvalid.pid = {pvalid.pid[np.where(res > 1)[0]]}\n" 
                                 f" pvalid.xdepth = {pvalid.xdepth[np.where(res > 1)[0]]}\n"
                                 f" pvalid.xdepth_inter = {pvalid.xdepth_inter[np.where(res > 1)[0]]}\n"
                                 f" pvalid.xdepth_decay = {pvalid.xdepth_decay[np.where(res > 1)[0]]}")
        
        self.final_stack.append(pvalid[np.where((pvalid.xdepth_inter >= max_xdepth) &
                                (pvalid.xdepth_decay >= max_xdepth))[0]])
        
        self.inter_stack = pvalid[np.where(pvalid.xdepth_inter < pvalid.xdepth_decay)[0]]
        self.decay_stack.append(pvalid[np.where(pvalid.xdepth_decay < pvalid.xdepth_inter)[0]])
        
        print(f"After set_xdepths:")
        print(f"working = {len(self.working_stack)}")
        print(f"interac = {len(self.inter_stack)}")
        print(f"decaying = {len(self.decay_stack)}")
        
        
        
    
    def run_interaction(self):
        self.children_stack.clear()
        self.failed_stack.clear()
        self.number_of_interactions += self.hadron_interaction.run_event_generator(
                                                    parents = self.inter_stack, 
                                                    children = self.children_stack, 
                                                    failed_parents = self.failed_stack)
        
        print("After run inter:")
        print(f"children = {len(self.children_stack)}")
        print(f"failded = {len(self.failed_stack)}")
        print(f"ninter = {self.number_of_interactions}")
        
              
        if len(self.failed_stack) > 0:
            fvalid = self.failed_stack.valid()
            finals_true = np.isin(fvalid.pid, self.pdg_lists.mceq_finals)
            finals_false = np.logical_not(finals_true)           
            self.final_stack.append(fvalid[np.where(finals_true)])
            decaying_particles = fvalid[np.where(finals_false)]
            decaying_particles.filter_code[:] = FilterCode.XD_DECAY_OFF.value
            self.decay_stack.append(decaying_particles)
            
            print("Distribute failed stack")
            print(f"append to finals = {len(finals_false)}")
            print(f"append to decay = {len(decaying_particles)}")
            
           
         
    def filter_particles(self):
        decay_pdgs = np.array([111], dtype = np.int32)
        
        cstack = self.children_stack.valid()
        
        decaying_true = np.isin(cstack.pid, decay_pdgs)
        self.decay_stack.append(cstack[np.where(decaying_true)[0]])
        
        finals_true = np.logical_or(np.isin(cstack.pid, self.pdg_lists.mceq_finals), 
                                    np.logical_and(cstack.energy < self.threshold_energy,
                                                   np.logical_not(np.isin(cstack.pid, decay_pdgs))))
        self.final_stack.append(cstack[np.where(finals_true)[0]])
        
        working_true = np.logical_not(np.logical_or(decaying_true, finals_true))
        self.working_stack.clear()
        print("After clear self.working_stack = ", len(self.working_stack))
        self.working_stack.append(cstack[np.where(working_true)[0]])
        print("After append self.working_stack = ", len(self.working_stack))
        
        print("From children:")
        print(f"to decaying_stack = {len(np.where(decaying_true)[0])}")
        print(f"to working_stack = {len(np.where(working_true)[0])}")
        print(f"to finals = {len(np.where(finals_true)[0])}")
        print(f"sum = {len(np.where(finals_true)[0]) + len(np.where(working_true)[0])+len(np.where(decaying_true)[0])}")
        print(f"cstack = {len(cstack)}")
        
        
    
    def decay_particles(self):
        if len(self.decay_stack) == 0:
            return
        self.working_stack.clear()
        check_xdepth(self.decay_stack)
        self.number_of_decays += self.decay_driver.run_decay(self.decay_stack, 
                                                             final_particles=self.working_stack,
                                                             stable_particles=self.final_stack)
        check_xdepth(self.final_stack)
        check_xdepth(self.working_stack)
        
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
    

        
        
    
    
    
    