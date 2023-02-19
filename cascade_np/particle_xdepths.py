from decay_xdepth import DecayXdepth
from csec_tables import CrossSectionTableMCEq, CSXdepthConversion
from pdg_pid_map import PdgLists
from csec_on_table import CrossSectionOnTable
from xdepth_on_table import XdepthOnTable
from xdepth_conversion import XdepthConversion
from MCEq.geometry.density_profiles import CorsikaAtmosphere
from tab_pproperties import TabulatedParticleProperties, ParticlePropertiesParticle

class NextDecayXdepth:
    def __init__(self, *, xdepth_conversion):   
        xdepth_on_table = XdepthOnTable(xdepth_conversion = xdepth_conversion, npoints=1000)
        particle_properties = ParticlePropertiesParticle()
        tab_particle_properties = TabulatedParticleProperties(particle_properties=particle_properties)
        self.decay_xdepth = DecayXdepth(tab_particle_properties=tab_particle_properties,
                                        xdepth_on_table=xdepth_on_table)
         
    def get_xdepth(self, pstack):        
        pslice = slice(0, len(pstack))
        pdg = pstack.pid[pslice]
        energy = pstack.energy[pslice]
        xdepth = pstack.xdepth[pslice]
        pstack.xdepth_decay[pslice] = (self.decay_xdepth
                                       .get_xdepth(pdg, energy, xdepth))
        
class NextInterXdepth:
    def __init__(self, *, xdepth_conversion):        
        cs_xdepth_conv = CSXdepthConversion()
        cs_table = CrossSectionTableMCEq(interaction_model="DPMJETIII191",
                                         cs_xdepth_conv = cs_xdepth_conv)
        cs_table.add_pdgs(PdgLists().longer_pi0_to_mceq)
        self.inter_xdepth = CrossSectionOnTable(cs_table)
        self.max_xdepth = xdepth_conversion.get_max_xdepth()
         
    def get_xdepth(self, pstack):       
        pslice = slice(0, len(pstack))
        pdg = pstack.pid[pslice]
        energy = pstack.energy[pslice]
        xdepth = pstack.xdepth[pslice]
        result_with_infs = (self.inter_xdepth
                                       .get_xdepth(pdg, energy) + xdepth)
                       
        pstack.xdepth_inter[pslice] = np.where(result_with_infs == np.inf, 
                                               self.max_xdepth, result_with_infs)
        
        
if __name__ == "__main__":
    from particle_array import ParticleArray
    import numpy as np
    
    atmosphere = CorsikaAtmosphere("SouthPole", "December")
    xconv =  XdepthConversion(atmosphere = atmosphere)
    xconv.set_theta(30)
    
    next_decay = NextDecayXdepth(xdepth_conversion=xconv)
    next_inter = NextInterXdepth(xdepth_conversion=xconv)
    
    pstack = ParticleArray(size=100)
    pstack.push(
        pid=np.array([111, 22, 111, 13, -13, 2212]),
        energy=np.array([2e10, 2e0, 2e10, 1e0, 1e0, 1e0]),
        xdepth=np.array([100, 56, 100, 98, 56, 500]),
    )
    
    print(pstack.xdepth_decay[0:len(pstack)])
    print(pstack.xdepth_inter[0:len(pstack)])
    
    next_decay.get_xdepth(pstack)
    next_inter.get_xdepth(pstack)
    
    print("After")
    print(pstack.xdepth_decay[0:len(pstack)])
    print(pstack.xdepth_inter[0:len(pstack)])
    
    
         