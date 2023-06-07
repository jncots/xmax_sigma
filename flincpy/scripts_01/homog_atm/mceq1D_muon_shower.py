import numpy as np
import mceq_config as config
from MCEq.core import MCEqRun
from MCEq.geometry.density_profiles import GeneralizedTarget
import matplotlib.pyplot as plt


config.debug_level = 1
config.max_density = 1.293e-3 #0.001225
config.e_min = 1e-1
config.muon_helicity_dependence = False
config.average_loss_operator = False
config.enable_2D = False
config.mceq_db_fname = ("/hetghome/antonpr/MCEq/MCEq/data/" + 
                        "mceq_db_lext_dpm191_v150.h5")


target_grammage = 10000 # g/cm2
target = GeneralizedTarget(len_target=target_grammage/config.max_density, 
                           env_density=config.max_density, 
                           env_name='homogeneous air')


mceq = MCEqRun(interaction_model="DPMJETIII191", 
               theta_deg=None, 
               density_model=target, 
               primary_model=None)


from normalized_spectrum import NormalizedSpectrum
from normalized_spectrum import HistFromDist

muon_spectr = NormalizedSpectrum(pdg_id = 13,
                                 etot_min = 0.5, etot_max = 1000, 
                                 spectral_index = 3,
                                 number_particle_norm = 1)

muon_init_spectrum = muon_spectr.dN_dekin(mceq.e_grid)
mceq.set_initial_spectrum(muon_init_spectrum, 
                          pdg_id=13, append=False)

mceq_xdepths = [0, 1000, 5000, 10000]
mceq.solve(int_grid=mceq_xdepths)


mceq_sol = {}
particle_groups = [(13, "mu+", "mu-"),
                    (14, "numu", "antinumu"),
                    (12, "nue", "antinue"),
                    (211, "pi+", "pi-"),
                    (11, "e+", "e-")]

for pdg, *group in particle_groups:
    
    dist_xdepth = {}
    for ixdepth in range(len(mceq_xdepths)):
        dist = None
        for particle in group:
            if dist is None:
                dist = mceq.get_solution(particle, grid_idx = ixdepth)
            else:
                dist += mceq.get_solution(particle, grid_idx = ixdepth)
        dist_xdepth[ixdepth] = dist
    
    mceq_sol[pdg] = dist_xdepth   
    
    
mceq_sol_hist = {}
for pdg in mceq_sol:
    xdepth_hist = {}
    for ixdepth in mceq_sol[pdg]:
        hist = HistFromDist(mceq.e_grid, mceq_sol[pdg][ixdepth])
        bins, dN = hist.hist(energy_bins)
        xdepth_hist[ixdepth] = dN
    mceq_sol_hist[pdg] = xdepth_hist    