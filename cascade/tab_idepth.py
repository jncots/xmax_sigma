import chromo
# from particle import Particle
import numpy as np
from tab_particles import TabParticles

tp = TabParticles()

# sibyll_valid_pdgs = ([113, 130, 211, 310, 321, 411, 421, 431,
#                       2112, 2212, 3112, 3122, 3212, 3222,
#                       3312, 3322, 4122, 4132, 4232, 4332])


# longer_pi0 = ([130, 211, 310, 321, 411, 421, 431,
#                511, 521, 531, 541, 2112, 2212, 3112,
#                3122, 3222, 3312, 3322, 3334, 4122,
#                4132, 4232, 4332, 5122, 5132, 5232, 5332])

# longer_30ps = ([130, 211, 310, 321, 2112, 2212, 3112, 3122, 3222, 3312, 3322, 3334])

# finals = ([11, 12, 13, 14, 16, 22])

# for i in tp._pdg2id:
# print(tp._id2pdg)

# event_kin = chromo.kinematics.FixedTarget(1e8, 2212, 2212)
# event_generator = chromo.models.DpmjetIII193(event_kin)

ist = 1

projectiles = []
# for i, pdg in enumerate(tp._id2pdg):

all_particles = tp._all_particles_dict

def get_projectiles():
    projectiles = []
    min_tau = 1e-19 # sec
    min_pdg = 22
    max_pdg = 6000
    
    
    for pdg in all_particles:        
        if all_particles[pdg].ctau is not None:
            tau = tp._all_particles_dict[pdg].ctau * 1e-1/2.998e10
        else:
            tau = np.nan
        
        if  (tau == np.inf or tau > min_tau) and (min_pdg < abs(int(pdg)) < max_pdg):
            projectiles.append(pdg)
        
    return projectiles    



# print(get_projectiles())

event_kin = chromo.kinematics.FixedTarget(1e8, 2212, 2212)
event_generator = chromo.models.Sibyll23d(event_kin)

projectiles = get_projectiles()
event_generator._projectiles = projectiles

# for pdg in get_projectiles():
#     print(pdg)
#     event_kin = chromo.kinematics.FixedTarget(1e2, int(pdg), 2212)
#     print(event_generator.cross_section(event_kin))
    
    

# def get_projectiles():

proj = []
for pdg in tp._all_particles_dict:
    name = tp._all_particles_dict[pdg].name
    
    if tp._all_particles_dict[pdg].ctau is not None:
        tau = tp._all_particles_dict[pdg].ctau * 1e-1/2.998e10
    else:
        tau = np.nan
        
    
    if tp._all_particles_dict[pdg].mass is not None:
        mass = tp._all_particles_dict[pdg].mass * 1e-3
    else:
        mass = np.nan        
    
    ctau = tau * 2.998e10 * 1e11/1e5
    # (1e-20 < tau < 1e-12)
    # if  (tau == np.inf or tau > 30e-12) and (22 < abs(int(pdg)) < 6000) :
    if  (tau == np.inf or tau > 30e-12) and (abs(int(pdg)) <= 22) :
        if mass is np.nan:
            mass = np.float64(0)  
        # print(f"{ist} {name}[{int(pdg)}], mass = {mass} GeV, tau = {tau} s, ctau = {ctau} km")
        print(f"{ist} {name}[{int(pdg)}],\t mass = {mass:.2e} GeV,\t tau = {tau:.2e} s, ctau = {ctau:.2e} km")
        ist += 1
        proj.append(int(pdg))
         
        
def get_sibyll_valid_pdgs():
    # Numbers from the code        
    valid_sib_ids = [7, 8, 9, 10, 11, 12, 13, 14, -13, -14]
    valid_sib_ids.extend([34, 35, 36, 37, 38, 39])
    valid_sib_ids.extend([59, 60, 71, 72, 74, 75])
    valid_sib_ids.extend([87, 88, 89, 99, 27])
    valid_pids = []
    for sib_id in valid_sib_ids:
        pdg_id = event_generator._lib.isib_pid2pdg(sib_id)
        valid_pids.append(abs(pdg_id))
        
    # valid_pids = list(set(valid_pids))
    valid_pids = list(set(valid_pids))
    valid_pids.sort()
       
    return valid_pids

vpdg = list(set([abs(p) for p in proj]))
vpdg.sort()
        
print(vpdg)

# for pdg in proj:
#     if abs(pdg) not in valid_pids:
#         print(pdg)
# # print(valid_pids)
# #             projectiles.append(pdg)      
    

# sybill_valid_pdgs
# interacting_pdgs
# final_pdgs
# sibyll_csection_map

    
    # try:
    #     event_kin = chromo.kinematics.FixedTarget(1e3, int(pdg), 2212)
    #     cs = event_generator.cross_section(event_kin)
    #     # print(event_generator.cross_section(event_kin))
    #     print(f"Valid pdg = {pdg}")
    # except Exception as ex:
    #     # pass
    #     print(ex)
    #     # print(f"Invalid pdg = {pdg}")    







# air = chromo.util.CompositeTarget([("N", 0.5), ("O", 0.5)])
# # (14, 7)

# def get_cross_section(target):
#     def total_cs(energy):
#         event_kin = chromo.kinematics.FixedTarget(energy, 2212, target)
#         try:
#             return np.float64(event_generator.cross_section(event_kin, precision=10000).inelastic)
#         except:
#             return np.float64(np.nan)
    
#     return total_cs 

# # pid_len = len(pdg_id_list)
# # energy_len = len(energy_array)

# sigma_tab = np.empty([len(pdg_id_list), len(energy_array)], dtype=np.float64)

# for pdg_id in pdg_id_list:
    
#     sigma_tab[_pdg2id[pdg], :] = np.frompyfunc(get_cross_section(pdg_id), 1, 1)(energy_array).astype("float64")

# sigma_tab[0, :] = np.array([20, 20])
# print(sigma_tab)


# event_kin = chromo.kinematics.FixedTarget(1e2, 2212, 2212)
# event_generator = chromo.models.Sibyll23d(event_kin)

# print(event_generator.cross_section())
# all_particles_list = Particle.findall()
# all_particles_dict = {p.pdgid: p for p in all_particles_list}

# plist = [int(p) for p in all_particles_dict.keys()]
# print([int(p) for p in all_particles_dict.keys()])


# event_kin = chromo.kinematics.FixedTarget(1e2, 2212, 2212)

# event_generator = chromo.models.Sibyll23d(event_kin)

# event_kin = chromo.kinematics.FixedTarget(1e2, 2212, 2212)

# print(event_generator.cross_section())

# event_kin = chromo.kinematics.FixedTarget(1e2, 2212, 2212)

# event_generator = chromo.models.Sibyll23d(event_kin)

# print(event_generator.cross_section())

# print([int(p) for p in all_particles_dict.keys()])


# def _sigma_wrapper(self, pid, energy):
#     event_kin = chromo.kinematics.FixedTarget(energy, pid, self.target)

#     if int(event_kin.p1) not in self.valid_pids:
#         print(f"{event_kin.p1} is not in valid_pids")
#         raise Exception("Not a valid pdg for beam")

#     self.event_generator.kinematics = event_kin
#     sigma_hair = self.event_generator._lib.sib_sigma_hair

#     pid_abs = abs(int(event_kin.p1))
#     if 1000 < pid_abs < 10000:
#         sigproj = 1
#     elif pid_abs % 1000 < 300:
#         sigproj = 2
#     else:
#         sigproj = 3
#     sigma = sigma_hair(sigproj, event_kin.ecm)
#     if isinstance(sigma, tuple):
#         return sigma[0]
#     return sigma




        