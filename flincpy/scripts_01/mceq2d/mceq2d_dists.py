from MCEq.core import MCEqRun
import mceq_config as config
import crflux.models as pm
from tqdm import tqdm

import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.integrate import nquad
from mceq2d_ihank import MCEqIHankel


def get_spline_2Dfunction(energy_grid, angle_grid, dist_function):
    # Create a spline representation of the input data
    spline = RectBivariateSpline(energy_grid, angle_grid, dist_function, kx = 5, ky = 5)
    
    # Define the 2D distribution function
    def distribution_function(energy, angle):
        return spline.ev(energy, angle)
    
    return distribution_function

def get_histogram(energy_bins, angle_bins, distribution_function):
    
    # distribution_function is assumed to be dN/(dE d(theta))
    
    # Initialize the histogram values to zero
    hist = np.zeros((len(energy_bins) - 1, len(angle_bins) - 1))

    # Integrate the distribution function in each bin
    for i in range(len(energy_bins) - 1):
        energy_start = energy_bins[i]
        energy_end = energy_bins[i+1]
        for j in range(len(angle_bins) - 1):
            angle_start = angle_bins[j]
            angle_end = angle_bins[j+1]
            hist[i,j], _ = nquad(distribution_function, 
                                 [(energy_start, energy_end), (angle_start, angle_end)])
            
    
    return energy_bins, angle_bins, hist


class MCEQDist2D():
    def __init__(self,
                 energy,
                 pdg_id,
                 theta_deg,
                 slant_depths,
                 #pname_tuples,
                 energy_range = [1e-1, 1e4],
                 interaction_model = "EPOS-LHC", 
                 hybrid_crossover = 0.1,
                 density_model = ("CORSIKA", ("USStd", None))):
        
        config.e_min = energy_range[0]
        config.e_max = energy_range[1]
        config.enable_2D = True
        config.mceq_db_fname = 'mceq_db_rare_decays_URQMD_lext_2D.h5'
        config.enable_default_tracking = False
        config.enable_em = False
        config.enable_em_ion = False
        config.hybrid_crossover = hybrid_crossover # 0.1
        config.muon_energy_loss = True
        config.enable_cont_rad_loss = True
        config.enable_energy_loss = True
        config.muon_helicity_dependence = True
        config.density_model = density_model

        config.adv_set['force_resonance'] = [421, 431, 411, 310]
        config.adv_set['disabled_particles'] = [22, 111, 16, 11]

        
        mceq_run = MCEqRun(
            #provide the string of the interaction model
            interaction_model=interaction_model,
            #primary cosmic ray flux model
            primary_model = (pm.HillasGaisser2012, "H3a"),
            # Zenith angle in degrees. 0=vertical, 90=horizontal
            theta_deg=theta_deg,
            density_model = density_model
        )

        self.mceq_run = mceq_run

        for slant_depth in slant_depths:
            if mceq_run.density_model.max_X < slant_depth:
                        raise ValueError(f"Maximum slant_xdepth = {mceq_run.density_model.max_X}")

        #Set the zenith angle
        mceq_run.set_theta_deg(theta_deg)
        mceq_run.set_single_primary_particle(energy, pdg_id = pdg_id)
        mceq_run.solve(int_grid=slant_depths)
        
        self.slant_depths = slant_depths
        
        self.spline_dist = None
        self.default_ebins = None
        self.default_angbins = None
        
    def set_particles(self, particles, idepth):
        # particles = [(-13, 0), (-13, -1), (-13, 1), (13, 0), (13, -1), (13, 1)]
        # part_tuple = [(-14, 0), (14, 0)]
    
        # Get results for each particle
        hank_trans_res = []
        for particle in particles:
            hank_trans_res.append(self.mceq_run.convert_to_theta_space(
                self.mceq_run.grid_sol, *particle))
            
        # Get a grid
        angle_grid = hank_trans_res[0][2]
        energy_grid = self.mceq_run.e_grid
        
        # Sum distributions for particles
        egrid_len = len(energy_grid)
        ang_dists = [None for i in range(egrid_len)]
        for i_energy in range(egrid_len):
            for res in hank_trans_res:
                if ang_dists[i_energy] is None:
                    ang_dists[i_energy] = res[3][0][i_energy]
                else:
                    ang_dists[i_energy] += res[3][0][i_energy]
        
        # Convert to numpy array           
        # ang_dists_dt is dN/dE (t*dt), where t is angle.
        ang_dists = np.array(ang_dists)
        # We multiply it to angle to get dN/dEdt        
        ang_dists_t = ang_dists*angle_grid[np.newaxis, :] 
        
        # Get spline of the function
        self.spline_dist = get_spline_2Dfunction(energy_grid, angle_grid, ang_dists_t)
        
        self.default_ebins = self.mceq_run.e_bins
        self.default_angbins = np.deg2rad(np.linspace(0, 90, 91))
        
        
    def histogram(self, energy_bins, angle_bins):   
        return get_histogram(energy_bins, angle_bins, self.spline_dist)
                

class MCEq2Histogram:
    def __init__(self, mceq_hankel, mceq_dist_2D, pdg):
        
        self.mceq_hankel = mceq_hankel
        self.mceq_run = mceq_dist_2D.mceq_run
        self.slant_depths = mceq_dist_2D.slant_depths
        
        if abs(pdg) == 13:
            particles = [(-13, 0), (-13, -1), (-13, 1), (13, 0), (13, -1), (13, 1)]
        else:
            particles = [(-pdg, 0), (pdg, 0)]  
        
        self.spline_dist = {}
        self.default_ebins = self.mceq_run.e_bins
        self.default_angbins = np.deg2rad(np.linspace(0, 90, 91))  
        
        self._set_particles(particles)
        
        
        
    def _set_particles(self, particles):
        # particles = [(-13, 0), (-13, -1), (-13, 1), (13, 0), (13, -1), (13, 1)]
        # part_tuple = [(-14, 0), (14, 0)]
    
        # Get results for each particle
        hank_trans_res = []
        
        (k_grid, 
        hankel_amp_part, 
        theta_grid, 
        theta_distr_part) = self.mceq_hankel.ihankel(particles)
        
                
        # Get a grid (for all particles and slant depths)
        angle_grid = theta_grid
        energy_grid = self.mceq_run.e_grid
        
        for idepth, slant_depth in enumerate(self.slant_depths):
            # Sum distributions for particles
            
            theta_distr_sum = None
            for particle_index in range(len(theta_distr_part)):
                if theta_distr_sum is None:
                    theta_distr_sum = theta_distr_part[particle_index][idepth, :, :]
                else:
                    theta_distr_sum += theta_distr_part[particle_index][idepth, :, :]
                    
                                
            # Convert to numpy array           
            # ang_dists_dt is dN/dE (t*dt), where t is angle.
            ang_dists = theta_distr_sum.T
            # We multiply it to angle to get dN/dEdt        
            ang_dists_t = ang_dists*angle_grid[np.newaxis, :] 
            
            # Get spline of the function
            self.spline_dist[slant_depth] = get_spline_2Dfunction(energy_grid, angle_grid, ang_dists_t)
        
       
    def histogram(self, energy_bins, angle_bins, slant_depth):   
        result =  get_histogram(energy_bins, angle_bins, self.spline_dist[slant_depth])
        hist_dict = {}
        hist_dict["en_bins"] = result[0]
        hist_dict["ang_bins"] = result[1]
        hist_dict["hist_en_ang"] = result[2]
        
        return hist_dict
    
    
class CalcMCEqHists:
    def __init__(self, mceq_sol, particles = [12, 13, 14]):
        self.particles = particles
        self.slant_depths = mceq_sol.slant_depths
        self.mceq_dists = {}
        
        #print("Calculating matrix for Hankel transformation...")
        mceq_hankel = MCEqIHankel(mceq_sol.mceq_run)
        #print("Finished")

        for particle in particles:
            self.mceq_dists[particle] = MCEq2Histogram(mceq_hankel, mceq_sol, particle)
            
            
        self.default_ebins = mceq_sol.mceq_run.e_bins
        self.default_angbins = np.deg2rad(np.linspace(0, 90, 91))       
    
    def histograms(self, energy_bins, angle_bins):
        mceq_hists = {}
        for particle in tqdm(self.particles, total=len(self.particles)):
            mceq_hist_depth = {}
            for slant_depth in self.slant_depths:
                print(f"Particle {particle}, slant depth = {slant_depth}")
                mceq_hist_depth[slant_depth] = (self.mceq_dists[particle].histogram(energy_bins, 
                                                                                    angle_bins,
                                                                                    slant_depth))
                
                
            mceq_hists[particle] = mceq_hist_depth
            
            
        return mceq_hists            
        
        
        
        
if __name__ == "__main__":
     
    mceq = MCEQDist2D(energy = 100,
                    pdg_id = 2212,
                    theta_deg = 30,
                    slant_depth = 143)
    
    
    