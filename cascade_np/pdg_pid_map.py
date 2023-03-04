import numpy as np


class PdgPidMap:
    def __init__(self, pid_pdg_dict, max_pdg = 6000):
        self.max_pdg = max_pdg
        self.pid_pdg_dict = pid_pdg_dict
        self._build_maps()
                         
    def _build_maps(self):
        max_pdg_in_dict = max([abs(pdg) for pdg in self.pid_pdg_dict])
        if  max_pdg_in_dict <= self.max_pdg:
            total_map = True
            max_pdg_map = max_pdg_in_dict
        else:
            total_map = False
            max_pdg_map = self.max_pdg
        
        pdg_pid = np.empty(len(self.pid_pdg_dict), dtype = np.int32)                   
        pid_pdg = np.full(2 * max_pdg_map + 1,-2147483640, dtype=np.int32)
        
        for pdg, pid in self.pid_pdg_dict.items():
            pdg_pid[pid] = pdg
            if abs(pdg) <= max_pdg_map:
                pid_pdg[pdg] = pid

        self.total_map = total_map
        self.max_pdg_map = max_pdg_map
        self.pdg_pid = pdg_pid
        self.pid_pdg = pid_pdg
        self.max_pid = max(pid_pdg)
            
    
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
            
            
        
class PdgLists:
    def __init__(self):
        self.sibyll_valid_pdgs = ([113, 130, 211, 310, 321, 411, 421, 431,
                        2112, 2212, 3112, 3122, 3212, 3222,
                        3312, 3322, 4122, 4132, 4232, 4332])

        self.longer_pi0 = ([130, 211, 310, 321, 411, 421, 431,
                511, 521, 531, 541, 2112, 2212, 3112,
                3122, 3222, 3312, 3322, 3334, 4122,
                4132, 4232, 4332, 5122, 5132, 5232, 5332])

        self.longer_30ps = ([130, 211, 310, 321, 2112, 2212, 
                        3112, 3122, 3222, 3312, 3322, 3334])

        self.final_pdgs = np.array([-11, 11, 12, -13, 13, 14, 16, 22], dtype = np.int32)
        self.finals = ([11, 12, 13, 14, 16, 22])
        

        self.mceq_hadrons = np.array([111, 130, -211, 211, 310,
                                       -321, 321, -2112, 2112, 
                                       -2212, 2212, -3122, 3122,
                                       ], dtype = np.int32)
            
        self.mceq_finals = np.array([-11, 11, 12, -13, 13, 14, 16, 22], dtype = np.int32)
        
        
        
        self.mceq_particles = np.array([-11, 11, 12, -13, 13, -14, 14, 16, 22,
                                        111, 130, -211, 211, 310,
                                       -321, 321, -2112, 2112, 
                                       -2212, 2212, -3122, 3122,
                                       ], dtype = np.int32)
        
        self.pid_pdg_mceq = ({111: 0, 130: 1, -211: 2, 211: 3, 
                    310: 4, -321: 5, 321: 6, -411: 7, 
                    411: 8, -421: 9, 421: 10, -431: 11, 
                    431: 12, -2112: 13, 2112: 14, -2212: 15, 
                    2212: 16, -3122: 17, 3122: 18})

        self._set_longer_pi0_to_mceq()
    
    
    def _set_longer_pi0_to_mceq(self):    
        self.longer_pi0_to_mceq = ({-511: -411, 511: 411, -521: -421, 521: 421,
                        -531: -431, 531: 431, -541: -431, 541: 431})
        
        for pdg in self.longer_pi0:
            if pdg > 3000:
                self.longer_pi0_to_mceq[-pdg] = -3122
                self.longer_pi0_to_mceq[pdg] = 3122
            
            