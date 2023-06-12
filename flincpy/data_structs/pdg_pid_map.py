import numpy as np
from particle import Particle

def unique_pdg_list(pdgs):
    res = list(set(pdgs))
    res.sort(key=lambda x: (abs(x), x > 0))
    return res

class PdgPidMap:
    def __init__(self, pid_pdg_dict, max_pdg = 6000):
        self.max_pdg = max_pdg
        self.pid_pdg_dict = pid_pdg_dict
        self._build_maps()
        
#     def _build_maps(self):
        
#         self.none_value = -2147483640
#         # Numpy array containing all known pdg ids
#         self.known_pdg_ids = np.array(unique_pdg_list(self.pid_pdg_dict.keys()))
        
#         max_pdg_in_dict = max([abs(pdg) for pdg in self.pid_pdg_dict])
#         if  max_pdg_in_dict <= self.max_pdg:
#             total_map = True
#             max_pdg_map = max_pdg_in_dict
#         else:
#             total_map = False
#             max_pdg_map = self.max_pdg
            
        
#         max_pid_in_dict = max([pid for pid in self.pid_pdg_dict.values()])   
             
        
#         pdg_pid = np.full(max_pid_in_dict + 2, self.none_value, dtype = np.int32)                   
#         pid_pdg = np.full(2 * max_pdg_map + 1, self.none_value, dtype=np.int32)
        
#         for pdg, pid in self.pid_pdg_dict.items():
#             pdg_pid[pid] = pdg
#             if abs(pdg) <= max_pdg_map:
#                 pid_pdg[pdg] = pid

#         self.total_map = total_map
#         self.max_pdg_map = max_pdg_map
#         self.pdg_pid = pdg_pid
#         self.pid_pdg = pid_pdg
#         self.max_pid = max(pid_pdg)
        
#         # Some element pointing to none_value
#         self.pdg_pid_default_ind =np.where(self.pdg_pid == self.none_value)[0][0]
#         self.pid_pdg_default_ind =np.where(self.pid_pdg == self.none_value)[0][0]    
                         
    def _build_maps(self):
        
        self.none_value = -2147483640
        # Numpy array containing all known pdg ids
        self.known_pdg_ids = np.array(unique_pdg_list(self.pid_pdg_dict.keys()))
        
        max_pdg_in_dict = max([abs(pdg) for pdg in self.pid_pdg_dict])
        if  max_pdg_in_dict <= self.max_pdg:
            total_map = True
            max_pdg_map = max_pdg_in_dict
        else:
            total_map = False
            max_pdg_map = self.max_pdg
        
        # Numpy array containing all known pdg ids
        self.known_pdg_ids = np.array(unique_pdg_list(self.pid_pdg_dict.keys()))
        
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
            
    
    def valid_pid_indices(self, pids):
        return np.where(np.isin(pids, self.known_pdg_ids))[0]
    
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
        

def pdg2mceq_idx_map(mceq_run):
    """Return dictionary with mapping from pdg to mceq_idx
    as mceq_run.pman.pdg2mceqidx gives

    Args:
        mceq_run (MCEq_Run): initialized MCEq_Run object
    """
    pdg_idx_map = {}
    for pdg_mceq, idx_mceq in mceq_run.pman.pdg2mceqidx.items():
        # Filter only "ordinary" particles
        if pdg_mceq[1] == 0 and abs(pdg_mceq[0]) < 10000 and idx_mceq > -1:  
            # print(pdg_mceq[0], idx_mceq)
            pdg_idx_map[pdg_mceq[0]] = idx_mceq
        
    return pdg_idx_map   
    
        
class PdgPidMap1:
    """Maps pdgs to pids and pids to pdgs
    """
    def __init__(self, pid_pdg_dict, max_pdg = 6000):
        self.max_pdg = max_pdg
        self.pid_pdg_dict = pid_pdg_dict
        self._build_maps()
                         
    def _build_maps(self):
        
        self.none_value = -2147483640
        
        max_pdg_in_dict = max([abs(pdg) for pdg in self.pid_pdg_dict])
        if  max_pdg_in_dict <= self.max_pdg:
            total_map = True
            max_pdg_map = max_pdg_in_dict
        else:
            total_map = False
            max_pdg_map = self.max_pdg
            
        
        max_pid_in_dict = max([pid for pid in self.pid_pdg_dict.values()])   
             
        
        pdg_pid = np.full(max_pid_in_dict + 2, self.none_value, dtype = np.int32)                   
        pid_pdg = np.full(2 * max_pdg_map + 1, self.none_value, dtype=np.int32)
        
        for pdg, pid in self.pid_pdg_dict.items():
            pdg_pid[pid] = pdg
            if abs(pdg) <= max_pdg_map:
                pid_pdg[pdg] = pid

        self.total_map = total_map
        self.max_pdg_map = max_pdg_map
        self.pdg_pid = pdg_pid
        self.pid_pdg = pid_pdg
        self.max_pid = max(pid_pdg)
        
        # Some element pointing to none_value
        self.pdg_pid_default_ind =np.where(self.pdg_pid == self.none_value)[0][0]
        self.pid_pdg_default_ind =np.where(self.pid_pdg == self.none_value)[0][0]
            
    
    def get_pdgs(self, pids):        
        try:
            return self.pdg_pid[pids]
        except:
            pids_np = np.array(pids)
            last_ind = len(self.pdg_pid) - 1
            valid_pids = np.where(abs(pids_np) < last_ind, 
                                  pids_np, 
                                  self.pdg_pid_default_ind)
            return self.pdg_pid[valid_pids]
                
    def get_pids(self, pdgs):
        try:
            return self.pid_pdg[pdgs]
        except:
            pdgs_np = np.array(pdgs)
            valid_pdgs = np.where(abs(pdgs_np) < self.max_pdg_map, 
                     pdgs_np, 
                     self.pid_pdg_default_ind)
            return self.pid_pdg[valid_pdgs]
            # return self._get_pids_full(pdgs)
    
    def _get_pids_full(self, pdgs):
        return np.array([self.pid_pdg_dict[pdg] for pdg in pdgs], 
                        dtype = np.int32)        
            
            
        
class PdgLists:
    def __init__(self):
        
        all_particles_dict = {p.pdgid: p for p in Particle.findall()}
        pdgs_6000 = [] 
        for p in all_particles_dict:
            pdg = int(p)
            if 10 < abs(pdg) < 6000:
                pdgs_6000.append(pdg)
                
        pdgs_6000.sort(key=lambda x: (abs(x), x > 0))
        self.pdgs_below_abs6000 = np.array(pdgs_6000, dtype = np.int32)
        
        
        
        
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
                
        self.mceq_particles = np.array([-11, 11, -12, 12, -13, 13, -14, 14, -16, 16, 22, 
                                        111, 130, -211, 211, 310, -321, 321, 
                                        -411, 411, -421, 421, -431, 431, 
                                        -2112, 2112, -2212, 2212, -3122, 3122
                                       ], dtype = np.int32)
        
        
        self.not_mceq_below_abs6000 = (self.pdgs_below_abs6000[np.where(np.logical_not(
            np.isin(self.pdgs_below_abs6000, self.mceq_particles)))[0]])      
        
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
            
            