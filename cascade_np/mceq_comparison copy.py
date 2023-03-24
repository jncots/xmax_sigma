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
                 pname_tuples,
                 interaction_model = "DPMJET-III-19.1", 
                 enable_energy_loss = True, 
                 muon_helicity_dependence = True,
                 disable_decays = []):
        
        config.mceq_db_fname = "mceq_db_lext_dpm191_v150.h5"
        config.enable_energy_loss = enable_energy_loss
        config.muon_helicity_dependence = muon_helicity_dependence
        config.adv_set["disable_decays"] = disable_decays
        
        print(json.dumps(config.adv_set, indent = 4))
        
        mceq_run = MCEqRun(
            #provide the string of the interaction model
            interaction_model=interaction_model,
            #primary cosmic ray flux model
            primary_model = (pm.HillasGaisser2012, "H3a"),
            # Zenith angle in degrees. 0=vertical, 90=horizontal
            theta_deg=theta_deg
        )

        #obtain energy grid (fixed) of the solution for the x-axis of the plots
        self.e_grid = mceq_run.e_grid
        self.e_widths = mceq_run.e_widths
        self.e_bins = mceq_run.e_bins

        #Set the zenith angle
        mceq_run.set_theta_deg(theta_deg)
        n_pts = 100
        X_grid = np.linspace(0.1, mceq_run.density_model.max_X, n_pts)
        mceq_run.set_single_primary_particle(energy, pdg_id = pdg_id)
        mceq_run.solve(int_grid=[slant_depth])

        # Populate longitudinal spectra for all particles:
        part_long_spectra = {}
        for p in mceq_run.pman.all_particles:
            print(f"particle = {p.name}")
            longitudinal_spectrum = []
            for idx in range(n_pts):
                longitudinal_spectrum.append(mceq_run.get_solution(p.name, grid_idx=0))
            
            part_long_spectra[p.name] = (p, longitudinal_spectrum)
    
    
        X_bins = X_grid + (X_grid[1]+X_grid[0])/2
        X_bins = np.array([0e0] + list(X_bins))
        xgrid_idx = np.digitize([slant_depth], X_bins)[0]
        
        print(f"Xbins = {X_bins[xgrid_idx - 1], X_bins[xgrid_idx]}")
        
        
        self.flux = dict()
        
        for pnames in pname_tuples:
            
            group_name = pnames[0]
            self.flux[group_name] = None
            for pname in pnames[1:]:
                if self.flux[group_name] is None:
                    self.flux[group_name] = part_long_spectra[pname][1][xgrid_idx]
                else:
                    self.flux[group_name] += part_long_spectra[pname][1][xgrid_idx] 
                
            self.flux[group_name] = self.flux[group_name] * self.e_widths

            
        
        # self.mu_spec =(e_grid, (part_long_spectra["mu+"][1][xgrid_inx[0]]+
        #             # part_long_spectra["mu+_l"][1][xgrid_inx[0]]+
        #             # part_long_spectra["mu+_r"][1][xgrid_inx[0]]+
        #             part_long_spectra["mu-"][1][xgrid_inx[0]]
        #             # part_long_spectra["mu-_r"][1][xgrid_inx[0]]+
        #             # part_long_spectra["mu-_l"][1][xgrid_inx[0]]
        #             )*e_width, r"${\mu}^{+} + {\mu}^{-}$ mceq")

        # self.numu_spec = (e_grid, (part_long_spectra["numu"][1][xgrid_inx[0]]+
        #             part_long_spectra["antinumu"][1][xgrid_inx[0]])*e_width, r"$\bar{\nu}_{\mu} + {\nu}_{\mu}$ mceq")
        
        
        # self.nue_spec = (e_grid, (part_long_spectra["nue"][1][xgrid_inx[0]]+
        #             part_long_spectra["antinue"][1][xgrid_inx[0]])*e_width, r"$\bar{\nu}_{e} + {\nu}_{e}$ mceq")

        # self.pi_spec = (e_grid, (part_long_spectra["pi+"][1][xgrid_inx[0]]+
        #             part_long_spectra["pi-"][1][xgrid_inx[0]])*e_width, r"$\bar{\pi}^{+} + {\pi}^{-}$ mceq")


        # self.ebins = mceq_run.e_bins
        # self.egrid = e_grid
    