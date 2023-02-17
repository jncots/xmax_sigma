import numpy as np
from pdg_array_map import PdgArrayMap

                          
class CrossSectionOnTable:
    def __init__(self, *, sigma_table, energy_grid, pid_pdg_dict):
        self.sigma_table = sigma_tab
        self.energy_grid = energy_grid
        self.pmap = PdgArrayMap(pid_pdg_dict)

    def get_pid_array(self, pdg_array):    
        return np.frompyfunc(self.pdg2pid, 1, 1)(pdg_array).astype("int32")
    
    def cross_section(self, pdg, energy, xdepth):
        pid_array = self.pmap.get_pids(pdg)
        for pid in set(pid_array):
            pid_slice = np.where(pid_array == pid)[0]
            print(pid_slice)
            xdepth[pid_slice] = np.interp(energy[pid_slice], self.energy_grid, self.sigma_tab[pid,:])

 
if __name__ == "__main__":

    from xmax_sigma.cascade.cs_table_mceq import TabulateCSMCEq
    from pdg_array_map import PdgArrayMap
    
    mceq_cs_tab = TabulateCSMCEq()
    sigma_tab = mceq_cs_tab.get_sigma()
    energy_grid = mceq_cs_tab.get_energy_grid()
    pid_from_pdg = mceq_cs_tab.get_pid_from_pdg()

    pdg_list = np.array([111, 111, 111, 
                        111, -211, 
                        111, 111, 2212, 111], dtype=np.int32)
    energy_list = np.array([5e3, 5e2, 4e6, 4e3, 2e5, 2e5, 4e3, 5e2, 2e5], dtype=np.float64) 
    xdepth = np.empty([len(pdg_list)], dtype=np.float64)    

    tcs = CrossSectionOnTable(sigma_tab = sigma_tab, energy_grid = energy_grid, pid_pdg_dict = pid_from_pdg) 
    tcs.cross_section(pdg_list, energy_list, xdepth)

    print(xdepth)