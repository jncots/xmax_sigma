from MCEq.data import InteractionCrossSections, HDF5Backend
from MCEq.core import MCEqRun
import crflux.models as pm
import numpy as np


class CrossSectionTableMCEq:
    def __init__(self, interaction_model="DPMJETIII191"):
        self.interaction_model = interaction_model
        self.mceq_run = MCEqRun(
        #provide the string of the interaction model
        interaction_model=self.interaction_model,
        #primary cosmic ray flux model
        primary_model = (pm.HillasGaisser2012, "H3a"),
        # Zenith angle in degrees. 0=vertical, 90=horizontal
        theta_deg=0.0
        )

        hdf5_backend = HDF5Backend()
        self.interaction_cs = InteractionCrossSections(hdf5_backend, self.interaction_model)
        self.set_energy_grid()

    def set_energy_grid(self, *, emin=1e0, emax=1e11, npoints=1000):
        pid_size = 0
        for p in self.mceq_run.pman.all_particles:
            if p.is_hadron:
                pid_size += 1
                
        self.energy_grid = np.geomspace(emin, emax, npoints, dtype='float64')
        self.sigma_tab = np.empty([pid_size, len(self.energy_grid)], dtype=np.float64)
        self._tabulate()
        
    def _tabulate(self):
        pid = 0
        self.pid_pdg = dict()  
        for p in self.mceq_run.pman.all_particles:
            if p.is_hadron:
                pdg = p.pdg_id[0]
                self.sigma_tab[pid, :] = np.interp(self.energy_grid, 
                                            self.interaction_cs.energy_grid.c + p.mass, 
                                            self.interaction_cs.get_cs(p.pdg_id[0], True))
                
                self.pid_pdg[pdg] = pid
                pid += 1
            
    def get_sigma(self):
        return self.sigma_tab
    
    def get_energy_grid(self):
        return self.energy_grid
    
    def get_pid_pdg_dict(self):
        return self.pid_pdg





# pid_size = 0
# for p in mceq_run.pman.all_particles:
#     if p.is_hadron:
#         pid_size += 1
        
# energy_grid = np.geomspace(1e0, 1e11, 1000, dtype='float64')
# sigma_tab = np.empty([pid_size, len(energy_grid)], dtype=np.float64)


# pid = 0
# pid_from_pdg = dict()  
# for p in mceq_run.pman.all_particles:
#     if p.is_hadron:
#         pdg = p.pdg_id[0]
#         sigma_tab[pid, :] = np.interp(energy_grid, 
#                                       interaction_cs.energy_grid.c + p.mass, 
#                                       interaction_cs.get_cs(p.pdg_id[0], True))
        
#         pid_from_pdg[pdg] = pid
#         pid += 1
     