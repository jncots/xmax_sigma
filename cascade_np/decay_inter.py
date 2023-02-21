import chromo
import numpy as np
from pathlib import Path
from particle_array import ParticleArray, FilterCode
from particle_xdepths import default_xdepth_getter
# from particle_event import CascadeParticle

chormo_path = Path(chromo.__file__).parent


# class Pythia8DecayAfterburner:
#     def __init__(self):
#         self._init_pythia()
#         self.number_of_decays = 0

#     def _init_pythia(self):
#         import importlib
#         from random import randint

#         lib = importlib.import_module(f"chromo.models._pythia8")
#         xml_path = chormo_path / "iamdata/Pythia8/xmldoc"
#         self.pythia = lib.Pythia(str(xml_path), False)
#         seed = randint(1, 10000000)
#         self.pythia.settings.resetAll()
#         self.pythia.readString("Random:setSeed = on")
#         self.pythia.readString(f"Random:seed = {seed}")
#         self.pythia.readString("Print:quiet = on")
#         self.pythia.readString("ProcessLevel:all = off")
#         self.pythia.readString("ParticleDecays:tau0Max = 1e100")
#         self.pythia.init()




class DecayInteraction:
    def __init__(self):
        self._init_pythia()
        self.number_of_decays = 0
        
    def _init_pythia(self):
        import importlib
        from random import randint

        lib = importlib.import_module(f"chromo.models._pythia8")
        xml_path = chormo_path / "iamdata/Pythia8/xmldoc"
        self.pythia = lib.Pythia(str(xml_path), False)
        seed = randint(1, 10000000)
        self.pythia.settings.resetAll()
        self.pythia.readString("Random:setSeed = on")
        self.pythia.readString(f"Random:seed = {seed}")
        self.pythia.readString("Print:quiet = on")
        self.pythia.readString("ProcessLevel:all = off")
        self.pythia.readString("ParticleDecays:tau0Max = 1e100")
        self.pythia.init()
        self.pythia.particleData.mayDecay(211, True)
        self.pythia.particleData.mayDecay(-211, True)
        self.pythia.particleData.mayDecay(13, True)
        self.pythia.particleData.mayDecay(-13, True)
        self.pythia.particleData.mayDecay(111, True)
    
    
    def calc_xdepth(self, pstack):
        slice_to_calc = np.where(pstack.valid().filter_code == FilterCode.XD_DECAY_OFF.value)[0]
        # print("slice_to_calc = ", pstack.valid().filter_code)
        # print("slice_to_calc = ", slice_to_calc)
        stack_to_calc = pstack[slice_to_calc]
        
        default_xdepth_getter.get_decay_xdepth(stack_to_calc)
        pstack[slice_to_calc] = stack_to_calc
    
    
    def generations(self, parents, zero_gen_len, pstack):
        
        pp = np.empty(len(parents) + 1, dtype = np.int32)
        pp[:-1] = parents - 1
        pp[-1] = -1
        
        
        pp1 = np.copy(pp)
        for i in range(100):    
            sl1 = np.where((pp1 < zero_gen_len)&(pp1 > -1))[0]
            
            pstack.xdepth[sl1] = pstack.xdepth_decay[pp[sl1]]
            pstack.generation_num[sl1] = pstack.generation_num[pp[sl1]] + 1
            pstack.filter_code[sl1] = FilterCode.XD_DECAY_OFF.value
            self.calc_xdepth(pstack)
            
            # print("ppp = ", pp[sl1], sl1)
            # print("pstack.xdepth", pstack.xdepth_decay[pp[sl1]])
            # print("pstack.xdepth", pstack.xdepth[sl1])
            # print("pstack.xdepth_decay", pstack.xdepth_decay[sl1])
            pp1 = pp[pp1]
            if not len(sl1) > 0:
                break

            
        
    def __call__(self, pstack):
        
        self.calc_xdepth(pstack)
        self.pythia.event.reset()

        for ip in range(len(pstack)):
            m0 = self.pythia.particleData.findParticle(pstack.pid[ip]).m0
            self.pythia.event.append(
                pstack.pid[ip],
                91,
                0,
                0,
                0,
                0,
                np.sqrt((pstack.energy[ip] + m0) * (pstack.energy[ip] - m0)),  # pz
                pstack.energy[ip],
                m0,
            )
            # self.pythia.particleData.mayDecay(pstack.pid[ip], True)

        # Decay it
        self.pythia.forceHadronLevel()
        
        
        self.res_stack = ParticleArray(1000)     
        self.res_stack.push(pid = self.pythia.event.pid(),
                       energy = self.pythia.event.en())
        
        self.res_stack.valid().xdepth_decay[0:len(pstack)] = pstack.valid().xdepth_decay
        self.res_stack.valid().generation_num[0:len(pstack)] = pstack.valid().generation_num
        self.res_stack.valid().filter_code[0:len(pstack)] = pstack.valid().filter_code
        
        
        # print(self.pythia.event.pid())
        self.generations(self.pythia.event.parents()[:,0], len(pstack), self.res_stack)
        self.res_stack = self.res_stack[np.where(self.pythia.event.status() == 1)]    
    
    def get_finals(self):
        return self.res_stack
        

if __name__ == "__main__":
    
    dint = DecayInteraction()
    
    pstack = ParticleArray(10)
    
    pstack.push(pid = [3312, 211, 111], energy = [1e3, 1e3, 1e3], 
                xdepth = [1e1, 1e1, 1e1])
    
    pstack.valid().filter_code[:] = FilterCode.XD_DECAY_OFF.value
    
    
    print(pstack.valid().filter_code)
    dint.calc_xdepth(pstack)
    print(pstack.valid().xdepth)
    print(pstack.valid().xdepth_decay)
    print(pstack.valid().filter_code)
    
       