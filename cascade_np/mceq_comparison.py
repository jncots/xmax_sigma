from scipy import integrate
import matplotlib.pyplot as plt
import numpy as np

#import solver related modules
from MCEq.core import MCEqRun
import mceq_config as config
#import primary model choices
import crflux.models as pm


class GetMCEQDist():
    def __init__(self):
        mceq_run = MCEqRun(
        #provide the string of the interaction model
        interaction_model='SIBYLL2.3c',
        #primary cosmic ray flux model
        primary_model = (pm.HillasGaisser2012, "H3a"),
        # Zenith angle in degrees. 0=vertical, 90=horizontal
        theta_deg=30.0
        )

        #Power of energy to scale the flux (the results will be returned as E**mag * flux)
        mag = 0

        #obtain energy grid (fixed) of the solution for the x-axis of the plots
        e_grid = mceq_run.e_grid
        mceq_run.e_widths

        #Dictionary for results
        flux = {}

        #Define a zenith angle, counted positively from vertical direction. Theta = 0. means vertical, theta = 90. horizontal
        theta = 30

        #Set the zenith angle
        mceq_run.set_theta_deg(theta)
        n_pts = 100
        X_grid = np.linspace(0.1, mceq_run.density_model.max_X, n_pts)
        mceq_run.set_single_primary_particle(1e3, pdg_id = 2212)
        mceq_run.solve(int_grid=X_grid)

        # Populate longitudinal spectra for all particles:
        part_long_spectra = {}
        for p in mceq_run.pman.all_particles:
            longitudinal_spectrum = []
            for idx in range(n_pts):
                longitudinal_spectrum.append(mceq_run.get_solution(p.name, grid_idx=idx))
            
            part_long_spectra[p.name] = (p, longitudinal_spectrum)
    
    
        xbin1 = X_grid + (X_grid[1]+X_grid[0])/2
        xbin1 = np.array([0e0] + list(xbin1))
        # print(X_grid[1]-X_grid[0])
        # print(xbin1, X_grid)

        xd = 587
        xgrid_inx = np.digitize([xd], xbin1)
        print(np.digitize([xd], xbin1))
        # numpy.digitize(x, bins, right=False)[source]
        
        self.mu_spec =(e_grid, (part_long_spectra["mu+"][1][xgrid_inx[0]]+
                    part_long_spectra["mu+_l"][1][xgrid_inx[0]]+
                    part_long_spectra["mu+_r"][1][xgrid_inx[0]]+
                    part_long_spectra["mu-"][1][xgrid_inx[0]]+
                    part_long_spectra["mu-_r"][1][xgrid_inx[0]]+
                    part_long_spectra["mu-_l"][1][xgrid_inx[0]])*e_grid, r"${\mu}^{+} + {\mu}^{-}$ mceq")

        self.numu_spec = (e_grid, (part_long_spectra["numu"][1][xgrid_inx[0]]+
                    part_long_spectra["antinumu"][1][xgrid_inx[0]])*e_grid, r"$\bar{\nu}_{\mu} + {\nu}_{\mu}$ mceq")
        
        
        self.nue_spec = (e_grid, (part_long_spectra["nue"][1][xgrid_inx[0]]+
                    part_long_spectra["antinue"][1][xgrid_inx[0]])*e_grid, r"$\bar{\nu}_{e} + {\nu}_{e}$ mceq")

        # plt.semilogx(e_grid, (part_long_spectra["mu+"][1][xgrid_inx[0]]+
        #             part_long_spectra["mu+_l"][1][xgrid_inx[0]]+
        #             part_long_spectra["mu+_r"][1][xgrid_inx[0]]+
        #             part_long_spectra["mu-"][1][xgrid_inx[0]]+
        #             part_long_spectra["mu-_r"][1][xgrid_inx[0]]+
        #             part_long_spectra["mu-_l"][1][xgrid_inx[0]])*e_grid, label = r"${\mu}^{+} + {\mu}^{-}$")

        # plt.semilogx(e_grid, (part_long_spectra["numu"][1][xgrid_inx[0]]+
        #             part_long_spectra["antinumu"][1][xgrid_inx[0]])*e_grid, label = r"$\bar{\nu}_{\mu} + {\nu}_{\mu}$")
        # plt.semilogx(e_grid, (part_long_spectra["nue"][1][xgrid_inx[0]]+
        #             part_long_spectra["antinue"][1][xgrid_inx[0]])*e_grid, label = r"$\bar{\nu}_{e} + {\nu}_{e}$")
        # plt.xlim(1e-1, 1e3)
        # plt.grid()
        # plt.grid(visible=True, which='minor', linestyle='--')
        # plt.legend()
        # # plt.ylim(1e-7, 1e3)    