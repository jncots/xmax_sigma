import numpy as np
from particle_array import ParticleArray

class CascadeDriver:
    def __init__(self):
        self.working_stack = ParticleArray(10000)
        self.final_stack = ParticleArray(10000)
        self.decay_stack = ParticleArray(10000)
    
    
    def run(self, *, pdg, energy, xdepth = 0):
        self.working_stack.push(pid = pdg, 
                         energy = energy, 
                         xdepth = xdepth)
        
        # while len(self.working_stack) > 0:
        #     self.cascade_event.get_event_particles(self.working_stack)
        
        # self.final_stack.append(self.working_stack)
        
    # def filter_below_threshold(self):
    #     threshold_energy = 1e3
        
    #     wstack = self.working_stack.valid()
    #     finals_slice = np.where(wstack.energy < threshold_energy & wstack.valid_code == 1)[0]
        
    #     self.decay_stack.append(wstack[finals_slice])
    #     wstack.valid_code[finals_slice] = 0
        
           
    
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
        finals_slice = np.where((np.in1d(wstack.pid, final_pdgs) 
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
    cas_driver.run(pdg = 111, energy = 1e3)
    
    cas_driver.working_stack.push(pid = 2212)
    cas_driver.working_stack.push(pid = 11, energy = 789, xdepth = 0)
    cas_driver.working_stack.push(pid = np.array([-11, 342, 13]), energy = np.array([20, 20, 20]))
    
    cas_driver.filter_final_particles()
    cas_driver.filter_decaying_particles()
    
    print(cas_driver.working_stack.valid().pid)
    print(cas_driver.get_final_particles().valid().pid)
    print(cas_driver.get_decaying_particles().valid().pid)
    print(cas_driver.working_stack.valid().pid)
    print(cas_driver.working_stack.valid().valid_code)
    
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
    

        
        
    
    
    
    