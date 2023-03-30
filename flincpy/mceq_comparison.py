from scipy import integrate
import matplotlib.pyplot as plt
import numpy as np
import json

#import solver related modules
from MCEq.core import MCEqRun
import mceq_config as config
#import primary model choices
import crflux.models as pm


class MCEQDistributions():
    def __init__(self,
                 energy,
                 pdg_id,
                 theta_deg,
                 slant_depth,
                 energy_range,
                 pname_tuples,
                 interaction_model = "DPMJET-III-19.1", 
                 generic_losses_all_charged = True,
                 enable_energy_loss = True, 
                 muon_helicity_dependence = True,
                 disable_decays = [],
                 hybrid_crossover = 0.5,
                 density_model = ("CORSIKA", ("USStd", None))):
        
        # config.mceq_db_fname = "mceq_db_lext_dpm191_v150.h5"
        config.generic_losses_all_charged = generic_losses_all_charged
        config.enable_energy_loss = enable_energy_loss
        config.muon_helicity_dependence = muon_helicity_dependence
        config.adv_set["disable_decays"] = disable_decays
        config.hybrid_crossover = hybrid_crossover
        
        config.e_min = energy_range[0]
        config.e_max = energy_range[1]
        # # config.enable_2D = True
        # # config.mceq_db_fname = 'mceq_db_rare_decays_URQMD_lext_2D.h5'
        # config.enable_default_tracking = False
        # config.enable_em = False
        # config.enable_em_ion = False
        # config.hybrid_crossover = 0.1
        # config.muon_energy_loss = True
        # config.enable_cont_rad_loss = True
        # config.enable_energy_loss = True
        # config.muon_helicity_dependence = True
        # config.density_model = ("CORSIKA", ("USStd", None))

        # config.adv_set['force_resonance'] = [421, 431, 411, 310]
        # config.adv_set['disabled_particles'] = [22, 111, 16, 11]
        
        # print(json.dumps(config.adv_set, indent = 4))
        
        mceq_run = MCEqRun(
            #provide the string of the interaction model
            interaction_model=interaction_model,
            #primary cosmic ray flux model
            primary_model = (pm.HillasGaisser2012, "H3a"),
            # Zenith angle in degrees. 0=vertical, 90=horizontal
            theta_deg=theta_deg,
            density_model = density_model
        )

        #obtain energy grid (fixed) of the solution for the x-axis of the plots
        self.e_grid = mceq_run.e_grid
        self.e_widths = mceq_run.e_widths
        self.e_bins = mceq_run.e_bins

        if mceq_run.density_model.max_X < slant_depth:
            raise ValueError(f"Maximum slant_xdepth = {mceq_run.density_model.max_X}")

        #Set the zenith angle
        mceq_run.set_theta_deg(theta_deg)
        mceq_run.set_single_primary_particle(energy, pdg_id = pdg_id)
        mceq_run.solve(int_grid=[slant_depth])

        # Populate longitudinal spectra for all particles:
        part_long_spectra = {}
        for p in mceq_run.pman.all_particles:
            # print(f"Spectrum for {p.name}")
            try: 
                part_long_spectra[p.name] = mceq_run.get_solution(p.name, grid_idx=0)
            except Exception as ex:
                pass
                # print(ex)    

        self.flux = dict()
        for pnames in pname_tuples:
            
            group_name = pnames[0]
            self.flux[group_name] = None
            for pname in pnames[1:]:
                if self.flux[group_name] is None:
                    self.flux[group_name] = part_long_spectra[pname]
                else:
                    self.flux[group_name] += part_long_spectra[pname]
                
            self.flux[group_name] = self.flux[group_name] * self.e_widths
    