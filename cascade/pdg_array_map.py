import numpy as np

class PdgArrayMap:
    def __init__(self, pid_pdg_dict, max_pdg = 6000):    
        max_pdg_in_dict = max([abs(pdg) for pdg in pid_pdg_dict])
        if  max_pdg_in_dict <= max_pdg:
            total_map = True
            max_pdg_map = max_pdg_in_dict
        else:
            total_map = False
            max_pdg_map = max_pdg
        
        pdg_pid = np.empty(len(pid_pdg_dict), dtype = np.int32)                   
        pid_pdg = np.full(2 * max_pdg_map + 1, np.nan, dtype=np.int32)
        
        for pdg, pid in pid_pdg_dict.items():
            pdg_pid[pid] = pdg
            pid_pdg[pdg] = pid

        self.total_map = total_map
        self.max_pdg_map = max_pdg_map
        self.pdg_pid = pdg_pid
        self.pid_pdg = pid_pdg
        self.pid_pdg_dict = pid_pdg_dict
        
        
    # def _build_arrays(self, )    
    
    def get_pdgs(self, pids):
        return self.pdg_pid[pids]
                
    def get_pids(self, pdgs):
        try:
            return self.pid_pdg[pdgs]
        except:
            return self._get_pids_full(pdgs)
    
    def _get_pids_full(self, pdgs):
        return np.array([self.pid_pdg_dict[pdg] for pdg in pdgs], 
                        dtype = np.int32)
            
            
        
        
        
    