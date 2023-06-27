from data_structs.pdg_pid_map import PdgPidMap1, pdg2mceq_idx_map
import numpy as np

class MceqGridCollector:
    def __init__(self, mceq_run, int_grid):
        self.mceq_run = mceq_run
        self.int_grid = int_grid
        self._init_pid_grid()
        self._init_sdepth_grid()
        self._init_energy_grid()
        self._init_shower_dist()
    
    
    def _init_pid_grid(self):
        """Initialization of the arrays related to
        particle id grid in mceq
        """
        self.pdg2idx_mapper = PdgPidMap1(pdg2mceq_idx_map(self.mceq_run))
        self.n_particles = self.mceq_run.pman.n_cparticles
            
    
    def _init_sdepth_grid(self):
        """Initialization of the arrays related to
        particle slant depth grid in mceq
        """
        self.mceq_run._calculate_integration_path(int_grid=self.int_grid, 
                                                  grid_var="X")
        
        sdepth_grid = np.cumsum(self.mceq_run.integration_path[1])
        
        # Bin edges
        sdepth_bins = (sdepth_grid[:-1] + sdepth_grid[1:])/2
        # Left limit of first bin
        sdepth_bins = np.insert(sdepth_bins, 0, 0)
        # Right limit of last bin        
        max_X = 2*sdepth_grid[-1] - sdepth_bins[-1]
        sdepth_bins = np.append(sdepth_bins, max_X)
        
        sdepth_widths = sdepth_bins[1:] - sdepth_bins[:-1]
        
        self.sdepth_grid = sdepth_grid
        self.sdepth_bins = sdepth_bins
        self.sdepth_widths = sdepth_widths
        # Size of sdepth grid
        self.n_sdepth_grid = len(self.sdepth_grid)
            
    def _init_energy_grid(self):
        """Initialization of the arrays related to
        particle energy grid in mceq
        """
        self.energy_bins = self.mceq_run.e_bins
        self.energy_grid = self.mceq_run.e_grid
        self.energy_widths = self.mceq_run.e_widths
        self.n_energy_grid = len(self.energy_grid)
        
        self._endist_matrix()        
    
    def _init_shower_dist(self):
        """Initialization of the array keeping
        result of cascade histograming
        """
        self.shower_dist = np.zeros(
               [self.n_sdepth_grid, 
                self.n_particles,
                self.n_energy_grid], 
                dtype=np.float64)
        
    
    def _sdepth_bins(self, sdepth_start = 0, sdepth_end = 1095):
        
        self.mceq_run._calculate_integration_path(int_grid=None, grid_var="X")
        
        sdepth_grid = np.cumsum(self.mceq_run.integration_path[1])
        sdepth_bins = (sdepth_grid[:-1] + sdepth_grid[1:])/2
        
        sdepth_bins = np.insert(sdepth_bins, 0, sdepth_start)
        sdepth_bins = np.append(sdepth_bins, sdepth_end)
        sdepth_widths = sdepth_bins[1:] - sdepth_bins[:-1]
        return sdepth_grid, sdepth_bins, sdepth_widths
    
    def _filter_valid(self, batch):
        mceq_idx = self.pdg2idx_mapper.get_pids(batch.pid)
        
        # Filter particles with pdgs present in mceq
        valid_mceq_batch = batch[np.where(mceq_idx != self.pdg2idx_mapper.none_value)[0]]
        
        # Filter particles with sdepth inside sdepth_bins
        sdepth_idx = np.digitize(valid_mceq_batch.xdepth, 
                                 self.sdepth_bins, right=True)

        valid_sdepth_batch_idx = np.where((sdepth_idx > 0) &
                                          (sdepth_idx < len(self.sdepth_bins)))[0]
            
        valid_sdepth_batch = valid_mceq_batch[valid_sdepth_batch_idx]
        
        # Filter particles with energy inside energy_bins
        energy_idx = np.digitize(valid_sdepth_batch.energy, 
                                 self.energy_bins, right=True)
        
        valid_energy_batch_idx = np.where((energy_idx > 0) &
                                          (energy_idx < len(self.energy_bins)))[0]
        
        valid_energy_batch = valid_sdepth_batch[valid_energy_batch_idx]
        
        return valid_energy_batch
    
    
    
    def _indices_for_addition(self, batch):
        # Calculate indices along sdepth dimension
        sdepth_idx = np.digitize(batch.xdepth, 
                                 self.sdepth_bins, right=True)
        # subtract 1 to get index out of bin number
        sdepth_idx = sdepth_idx - 1
        
    
        # Calculate indices along energy-particle dimension
        energy_idx = np.digitize(batch.energy, 
                                 self.energy_bins, right=True)

        # subtract 1 to get index out of bin number
        energy_idx = energy_idx - 1
        mceq_idx = self.pdg2idx_mapper.get_pids(batch.pid)
        
        self.sdepth_idx = sdepth_idx
        self.energy_idx = energy_idx
        self.mceq_idx = mceq_idx
        self.batch_energy = batch.energy
        
        
    def _endist_matrix(self):
        """calculates `len(self.energy_grid)` matrices with 3x3 shape
        to find an energy distribution on `energy_grid`
        for a single particle. 
        
        The energy distribution distribute the particle in 3 bins
        on a grid so that it saves number, energy and energy**2 when
        dn/dE*self.energy_widths
        """
        emats = []
        for i in range(1, len(self.energy_grid)-1):
            arr_rows = []
            for j in range(3):
                arr_rows.append(self.energy_grid[i-1:i+2]**j)
            
            emats.append(np.linalg.inv(np.vstack(arr_rows) 
                                       * self.energy_widths[i-1:i+2]))
            
        emats.insert(0, emats[0])  
        emats.insert(-1, emats[-1])
        self.emats = np.array(emats)
        
        # Shift vector
        shift_idx = np.zeros(len(self.energy_grid), dtype = np.int32)
        shift_idx[1:-1] = -1
        shift_idx[-1] = -2
        self.shift_idx = shift_idx
        
    def add_particles(self, batch, nruns=None):
        """Add batch of particles to the grid (sdepth, pid, energy)

        Args:
            batch (ParticleArray): particles to add
            nruns (int, optional): number of runs (cascades)
            for the batch
        """
        valid_batch = self._filter_valid(batch)
        self._indices_for_addition(valid_batch)
        
        batch_dist = np.zeros(
               [self.n_sdepth_grid, 
                self.n_particles,
                self.n_energy_grid], 
                dtype=np.float64)
        
        for i in range(len(self.energy_grid)):
            pidx = np.where(self.energy_idx == i)[0]
            en = self.batch_energy[pidx]
            evec = np.vstack([en**0, en, en**2])
            weights = np.matmul(self.emats[i], evec)
            
            for j in range(3):
                np.add.at(batch_dist, 
                          (self.sdepth_idx[pidx], 
                            self.mceq_idx[pidx], 
                            self.energy_idx[pidx] + j + self.shift_idx[i]), 
                    weights[j, :])
            
        if isinstance(nruns, int) and nruns > 0:
            self.shower_dist += batch_dist/nruns
        else:
            self.shower_dist += batch_dist
                      
    def state_vectors(self):
        shape = self.shower_dist.shape
        return self.shower_dist.reshape(shape[0], shape[1]*shape[2])
                
        
    def shower_on_grid(self, pdg = None):
        """Return the resulting distribution on
        (sdepth, p_id, energy)

        Args:
            pdg (list): list to filter only for required pdgs
        """
        if pdg is None:
            return self.shower_dist
        else:
            # Return only part of a histogram 
            # for pdgs in `pdg` list
            pdg_np = np.unique(pdg)
            mceq_idx = self.pdg2idx_mapper.get_pids(pdg_np)
            
            valid_pdg_idx = np.where(mceq_idx != self.pdg2idx_mapper.none_value)[0]
            valid_pdg = pdg_np[valid_pdg_idx]
            valid_mceq_idx = mceq_idx[valid_pdg_idx]
                      
            return (self.shower_dist[:, valid_mceq_idx, :], # dist
                    valid_pdg,                   # list of valid pdgs
                    valid_mceq_idx)              # list of valid mceq particle indicies
        

        
        