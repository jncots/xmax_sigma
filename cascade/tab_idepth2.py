import chromo
# from particle import Particle
import numpy as np
from tab_particles import TabParticles

tp = TabParticles()

sibyll_valid_pdgs = ([113, 130, 211, 310, 321, 411, 421, 431,
                      2112, 2212, 3112, 3122, 3212, 3222,
                      3312, 3322, 4122, 4132, 4232, 4332])


longer_pi0 = ([130, 211, 310, 321, 411, 421, 431,
               511, 521, 531, 541, 2112, 2212, 3112,
               3122, 3222, 3312, 3322, 3334, 4122,
               4132, 4232, 4332, 5122, 5132, 5232, 5332])

longer_30ps = ([130, 211, 310, 321, 2112, 2212, 3112, 3122, 3222, 3312, 3322, 3334])

finals = ([11, 12, 13, 14, 16, 22])


# p[2212]
# n[2112]
# pi+[211]
# K+[321]


sibyll_cs_channel = [2112, 211, 321]
sibyll_cs_map = dict()
for pdg in sibyll_valid_pdgs:
    apdg = abs(pdg)
    if apdg == 2212 or apdg == 2112:
        sibyll_cs_map[pdg] = 2212
    elif apdg == 211:
        sibyll_cs_map[pdg] = 211
    else:
        sibyll_cs_map[pdg] = 321
                
        


event_kin = chromo.kinematics.FixedTarget(1e8, 2212, 2212)
event_generator = chromo.models.Sibyll23d(event_kin)

air = chromo.util.CompositeTarget([("N", 0.78), ("O", 0.22)])


def get_cross_section(proj_pdg, target):
    def total_cs(energy):
        event_kin = chromo.kinematics.FixedTarget(energy, proj_pdg, target)
        try:
            return np.float64(event_generator.cross_section(event_kin).total)
        except:
            return np.float64(np.nan)
    return total_cs 

energy_array = np.geomspace(1e0, 1e10, 10, dtype='float64')
sigma_tab = np.empty([len(sibyll_cs_channel), len(energy_array)], dtype=np.float64)

for i, pdg in enumerate(sibyll_cs_channel):
    sigma_tab[i, :] = np.frompyfunc(get_cross_section(pdg, 2212), 1, 1)(energy_array).astype("float64")
    
    
    
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




        