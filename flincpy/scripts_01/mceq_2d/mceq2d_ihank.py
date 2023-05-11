import mceq_config as config
from pyhank import HankelTransform

import numpy as np
from scipy.interpolate import interp1d, splrep


class MCEqIHankel:
    """The class provides a way to perform an inverse Hankel transform on a set of 
    particle distributions obtained from the MCEq2D code. The transform is performed 
    for a set of particles specified by their particle ID and helicity.

    The class has the following methods:

    init(self, mceq_run, k_grid_size = 1024):
    Initializes the class by setting the MCEq run object and the size of 
    the k-grid for the Hankel transform.

    set_ksize(self, k_grid_size):
    Sets the size of the k-grid for the Hankel transform.

    set_mceq(self, mceq_run):
    Sets the MCEq run object.

    _pick_particles(self, pgd_hels):
    Picks the solutions for the particles specified by their particle ID and helicity.

    _to_theta_grid(self, theta_grid, inverse_hankel_transfs):
    Interpolates the inverse Hankel transform to a given theta grid.

    _subdivide_for_particles(self, pgd_hels, theta_distr, hankel_amps):
    Divides the matrices for each particle.

    ihankel(self, pdg_hels, theta_grid = None):
    Performs the inverse Hankel transform for the particles 
    specified by their particle ID and helicity. The results are returned on a given theta grid or on the internal theta grid returned by the inverse Hankel transform. The method returns a tuple containing the k-grid, the Hankel amplitudes, the theta grid, and the particle distributions in theta space.
    
    Example:
    ### Initialize MCEqIHankel object
    mceq_hankel = MCEqIHankel(mceq_run)

    ### Perform inverse Hankel transform for muon neutrinos and antineutrinos with helicity 0
    k_grid, hankel_amp_part, theta_grid, theta_distr_part = mceq_hankel.ihankel([(14, 0), (-14, 0)])
    
    ### Increase k-grid size (increase accuracy, but reduce speed)
    mceq_hankel.set_ksize(2048)
    
    ### Change the mceq_run object:
    mceq_hankel.set_mceq(mceq_run)
    
    """
    def __init__(self, mceq_run, k_grid_size = 2048):
        self.set_mceq(mceq_run)
        self.set_ksize(k_grid_size)
        
    def set_ksize(self, k_grid_size):
        self.k_grid = np.linspace(np.min(config.k_grid), np.max(config.k_grid), k_grid_size)
        self.ht_obj = HankelTransform(order=0, k_grid = self.k_grid)
        
    def set_mceq(self, mceq_run):
        self.mceq_run = mceq_run    
    
    def _pick_particles(self, pgd_hels):
        len_egrid = len(self.mceq_run.e_grid)
        en_inds = []
        for pgd_hel in pgd_hels:
            start_ind =  self.mceq_run.pman.pdg2mceqidx[pgd_hel] * len_egrid
            end_ind = start_ind + len_egrid
            en_inds.append(np.arange(start_ind, end_ind))

        en_slice = np.concatenate(en_inds)
        return self.mceq_run.grid_sol[:,:,en_slice]
    
    def _to_theta_grid(self, theta_grid, inverse_hankel_transfs):
        
        if theta_grid is None:
            theta_grid = self.ht_obj.r
            return theta_grid, inverse_hankel_transfs
        else:
            theta_distr = interp1d(self.ht_obj.r, 
                                     inverse_hankel_transfs, 
                                     axis = 1,
                                     fill_value="extrapolate",
                                     kind = "cubic",
                                     )(theta_grid)
            return theta_grid, theta_distr
        
    def _subdivide_for_particles(self, pgd_hels, theta_distr, hankel_amps):
        len_egrid = len(self.mceq_run.e_grid)
        start_ind = 0
        theta_distr_part = []
        hankel_amps_part = []
        for i in enumerate(pgd_hels):
            end_ind = start_ind + len_egrid
            theta_distr_part.append(theta_distr[:,:, start_ind:end_ind])
            hankel_amps_part.append(hankel_amps[:,:, start_ind:end_ind])
            start_ind = end_ind
            
        return hankel_amps_part, theta_distr_part
    
    def ihankel(self, pdg_hels, theta_grid = None):
        """Make a inverse hankel transformation (to theta space) for particles
        in pdg_hels list on a theta_grid. If theta_grid is not given
        the results are returned on internal theta grid returned by inverse hankel transformation

        Args:
            pdg_hels (list): e.g. [(13, 0), (13, -1), (14, 0)]
            theta_grid (ndarray): e.g. np.deg2rad(np.array(0, 90, 91))

        Returns:
            tuple: k_grid, hankel_amp_part, theta_grid, theta_distr_part
            
            theta_distr_part has the following structure:
                theta_distr_part[particle_index][slant_depth_index, theta_grid_index, mceq_energy_grid_index]
                
            hankel_amp_part has the following structure:
                hankel_amp_part[particle_index][slant_depth_index, k_grid_index, mceq_energy_grid_index]
        """
        
        # Pick solutions only for particles in pdg_hels list, 
        # e.g. [(13, 0), (13, -1), (14, 0)]
        hankel_transf = self._pick_particles(pdg_hels)
        
        # Resample hankel amplitudes to k grid of HankelTransform object
        k_grid = self.ht_obj.kr
        oversampled_hankel_amps = interp1d(config.k_grid, 
                                           hankel_transf, 
                                           axis = 1, 
                                           kind = "cubic",
                                           )(k_grid)
        
        # Note! kind = "linear" gives a larger flux at small angles
        # than kind = "cubic". It should be investigated further.
        # Number points in k_grid influence flux at large angles
        # Smaller number of points gives lower flux at large angles
        
    
        # Get inverse transformation (theta space) 
        # It should be multiplied to 2pi (not sure why, but it give the correct solution)
        inverse_hankel_transfs = self.ht_obj.iqdht(oversampled_hankel_amps, axis = 1)*(2*np.pi)
        
        # Interpolate to theta_grid
        theta_grid, theta_distr = self._to_theta_grid(theta_grid, inverse_hankel_transfs)
        
        # Divide the matrix to submatricies for each particle
        hankel_amp_part, theta_distr_part = self._subdivide_for_particles(pdg_hels, 
                                                         theta_distr, 
                                                         oversampled_hankel_amps)
        
        return (k_grid, 
                hankel_amp_part,
                theta_grid, 
                theta_distr_part
                )