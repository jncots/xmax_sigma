import chromo
import numpy as np
from pathlib import Path
from particle_array import ParticleArray, FilterCode

chormo_path = Path(chromo.__file__).parent

class DecayDriver:
    def __init__(self, xdepth_getter, decaying_pdgs = None, stable_pdgs = None):
        self.number_of_decays = 0
        self.xdepth_getter = xdepth_getter
        self.decaying_pdgs = decaying_pdgs
        self.stable_pdgs = stable_pdgs
        self._init_pythia()
        
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
        
        if self.stable_pdgs is not None:
            for pdg in self.stable_pdgs:
                self.pythia.particleData.mayDecay(pdg, False)
                
        if self.decaying_pdgs is not None:
            for pdg in self.decaying_pdgs:
                self.pythia.particleData.mayDecay(pdg, True)
    
    
    def set_xdepth_decay(self, pstack):
        """Fill xdepth_decay array if filter_code == FilterCode.XD_DECAY_OFF.value

        Args:
            pstack (ParticleArray): array of particle which need to set xdepth_decay
        """
        slice_to_fill = np.where(pstack.valid().filter_code == FilterCode.XD_DECAY_OFF.value)[0]
        # stack_to_fill is a copy, because of advanced indexing in numpy
        stack_to_fill = pstack[slice_to_fill]
        self.xdepth_getter.get_decay_xdepth(stack_to_fill)
        pstack[slice_to_fill] = stack_to_fill
    
    
    def fill_xdepth_for_decay_chain(self, pstack, parents, zero_generation_length):
        """The function follows the decay chain using information in "parent" array
        and fills `xdepth`, `xdepth_decay`, and `generation_num` of `pstack`
        It is assumed that 0th generation particles has already defined `xdepth_decay` 
        array
        
        Args:
            `pstack` (ParticleArray): array to fill
            `parents` (np.array): parent information
            `zero_generation_length` (int): length of 0th generation particles (the ones that decayed)
        """
        
        # parent_indices contains 0-based indices in pstack arrays
        # The last element is introduced for the indicies == -1
        # for "no parent" case
        parent_indices = np.empty(len(parents) + 1, dtype = np.int32)
        parent_indices[:-1] = parents - 1
        parent_indices[-1] = -1
        
        
        parent_gen = np.copy(parent_indices)
        # It is assumed that no more than 100 generations in the decay chain
        for _ in range(100):
            # Take elements with indicies smaller than the length of 0th generation
            # and filter out elements which point to "no parent"    
            generation_slice = np.where((parent_gen < zero_generation_length) & (parent_gen > -1))[0]
            # generation_slice contains elements for current generation (starting with 1st generation)
            # parent_indices[generation_slice] are corresponding indicies of parents
            pstack.xdepth[generation_slice] = pstack.xdepth_decay[parent_indices[generation_slice]]
            pstack.generation_num[generation_slice] = pstack.generation_num[parent_indices[generation_slice]] + 1
            # Set filter code to fill it in "set_xdepth_code()""
            pstack.filter_code[generation_slice] = FilterCode.XD_DECAY_OFF.value
            self.set_xdepth_decay(pstack)
            # parent_gen points to parents of parents ...
            parent_gen = parent_indices[parent_gen]
            # break of loop when next generation_slice is empty
            if not len(generation_slice) > 0:
                break

            
        
    def run_decay(self, pstack):
        """Run decay of particle in pstack
                
        FilterCode.XD_DECAY_OFF.value for `filter_code` should be set for particles 
        for which xdepth_decay is not set 


        Args:
            pstack (ParticleArray): stack with decaying particles
        """
        # Set xdepth_decay for particles which doesn't have it
        self.set_xdepth_decay(pstack)
        
        # Fill the Pythia stack of particles that should decay
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

        # Decay it
        self.pythia.forceHadronLevel()
        
        self.number_of_decays = len(np.where(self.pythia.event.status() == 2)[0])
        # Process event from Pythia
        decay_stack = ParticleArray(1000)     
        decay_stack.push(pid = self.pythia.event.pid(),
                       energy = self.pythia.event.en())
        
        # Set 0th generation
        gen0_slice = slice(0, len(pstack))
        decay_stack.valid().xdepth_decay[gen0_slice] = pstack.valid().xdepth_decay
        decay_stack.valid().generation_num[gen0_slice] = pstack.valid().generation_num
        decay_stack.valid().filter_code[gen0_slice] = pstack.valid().filter_code
        
        # Get parents array and fill in rest generations
        parents = self.pythia.event.parents()[:,0]        
        self.fill_xdepth_for_decay_chain(decay_stack, parents, len(pstack))
        # Filter final particles
        self.final_in_decay = decay_stack[np.where(self.pythia.event.status() == 1)]    
    
    def get_finals(self):
        return self.final_in_decay
        

if __name__ == "__main__":
    
    # Example:
    from particle_xdepths import default_xdepth_getter
    from particle_array import ParticleArray, FilterCode
    import numpy as np
    
    # Set array of particles to decay
    pstack = ParticleArray(10)    
    pstack.push(pid = [3312, 5232, 3312, 5232], energy = [1e3, 1e3, 1e2, 6e2], 
                xdepth = [1e1, 1e1, 1e1, 1e1])
    pstack.valid().filter_code[:] = FilterCode.XD_DECAY_OFF.value
    pstack.valid().generation_num[:] = np.int32(0)
    
    decay_driver = DecayDriver(default_xdepth_getter, 
                            decaying_pdgs=[111], 
                            stable_pdgs=[-211, 211, -13, 13],
                            # decaying_pdgs=[111, -211, 211, -13, 13],
                            )
    decay_driver.run_decay(pstack)
    print("pid = ", decay_driver.get_finals().valid().pid)
    print("energy = ", decay_driver.get_finals().valid().energy)
    print("xdepth_decay = ", decay_driver.get_finals().valid().xdepth_decay)
    print("xdepth = ", decay_driver.get_finals().valid().xdepth)
    print("gen_num = ", decay_driver.get_finals().valid().generation_num)
    print("number of finals = ", len(decay_driver.get_finals()))
    print("number of decays = ", decay_driver.number_of_decays)
    
       