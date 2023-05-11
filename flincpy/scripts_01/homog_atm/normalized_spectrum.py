import numpy as np
from particle import Particle


class NormalizedSpectrum:
    def __init__(self,
                 pdg_id,
                 etot_min,
                 etot_max, 
                 spectral_index = 3, 
                 number_particle_norm = 1):
        self.spectral_index = spectral_index
        self.number_particle_norm = number_particle_norm
        self.mass = Particle.from_pdgid(pdg_id).mass/1e3
        
        alpha = self.spectral_index - 1
        self.sp_norm = (self.number_particle_norm 
                * alpha/(1/etot_min**alpha - 1/etot_max**alpha))
        
        
        self.etot_min = etot_min
        self.etot_max = etot_max
        self._calc_limits()
    
    def _etot2ptot(self, etot, mass):
        return np.sqrt((etot - mass)*(etot + mass))
        
    def _calc_limits(self):
        self.ptot_min = self._etot2ptot(self.etot_min, self.mass)
        self.ptot_max = self._etot2ptot(self.etot_max, self.mass)
        
        self.ekin_min = self.etot_min - self.mass
        self.ekin_max = self.etot_max - self.mass
            
    
        
    def dN_dptot(self, ptot_grid): 
        result =  (self.sp_norm * ptot_grid/
                (ptot_grid**2 + self.mass**2)**((self.spectral_index + 1)/2)
                )
        
        return np.where((ptot_grid < self.ptot_min) | (ptot_grid > self.ptot_max), 
                                                       0, result)
        
        
    def dN_dekin(self, ekin_grid): 
        result =  self.sp_norm *(ekin_grid + self.mass)**(-self.spectral_index)
        return np.where((ekin_grid < self.ekin_min) | (ekin_grid > self.ekin_max), 
                                                       0, result)
        
    def dN_detot(self, etot_grid):
        result =  self.sp_norm *etot_grid**(-self.spectral_index)
        return np.where((etot_grid < self.etot_min) | (etot_grid > self.etot_max), 
                                                       0, result)
        
        
class HistFromDist:
    def __init__(self, xgrid, dNdx):
        self.xgrid = xgrid
        self.dNdx = dNdx
        
    def hist(self, xbins):
        
        xgrid_centers = (xbins[:-1] + xbins[1:])/2
        xgrid_widths = xbins[1:] - xbins[:-1]
        
        dNdx_centers = np.interp(x = xgrid_centers, 
                  xp = self.xgrid, 
                  fp = self.dNdx,
                  left = 0,
                  right= 0)
        
        dN = dNdx_centers * xgrid_widths
        
        return xbins, dN
     
    def normalized_hist(self, xbins):
        xbins, dN = self.hist(xbins)
        norm = 1.0/np.sum(dN)
        if norm == np.inf:
            norm = 0
        print(f"norm = {norm}")
        
        return xbins, norm*dN
        
        
            
        
                
        
        
                
           
        