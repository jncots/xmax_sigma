import impy
import random
import numpy as np
import dataclasses
import math

from abc import ABC, abstractmethod


@dataclasses.dataclass
class CascadeParticle:
    pid: int
    energy: float
    xdepth: float


class ParticleEvent(ABC):
    
    def __init__(self, particle):
        self.set_particle(particle)

    @abstractmethod
    def set_particle(self, particle):
        pass
        
    @abstractmethod    
    def get_xdepth(self):
        pass
        
    @abstractmethod    
    def get_products(self):
        pass


class HadronEvent(ParticleEvent):
    
    mbarn_in_cm2 = 1e-27
    proton_mass_g = 1.672621e-24
    mass_barn = proton_mass_g / mbarn_in_cm2
    
    
    def __init__(self, particle):
        super().__init__(particle)
        
        ekin = impy.kinematics.FixedTarget(10, "proton", (14, 7))
        self.event_generator = impy.models.Sibyll23d(ekin)

        self.set_average_A(14)
        self.valid_pids = self._get_valid_pids()
        
        
    def set_particle(self, particle):
        self.particle = particle
         
    def get_xdepth(self):
        """Returns slant depth of next interaction or None
        """
        
        self.last_xdepth = None
        try:
            average_xdepth = self.get_average_xdepth()
        except Exception:
            return None

        if not average_xdepth or math.isnan(average_xdepth):
            return None
          
        random_value = -np.log(1 - random.random())
        self.last_xdepth = random_value * average_xdepth
        return self.last_xdepth
           
    def get_products(self):       
        if not self.last_xdepth:
            return None
               
        total_xdepth = self.particle.xdepth + self.last_xdepth
        event = next(self.event_generator(1)).final_state()

        products = [] 
        for particle in event:
            products.append(CascadeParticle(particle.pid, particle.en, total_xdepth))    

        return products
    
    def _get_valid_pids(self):
        """Returns a set of valid pdg ids accepted by `Sibyll`
        event generator
        """
        valid_sib_ids = [7, 8, 9, 10, 11, 12, 13, 14, -13, -14]
        valid_sib_ids.extend([34, 35, 36, 37, 38, 39])
        valid_sib_ids.extend([59, 60, 71, 72, 74, 75])
        valid_sib_ids.extend([87, 88, 89, 99, 27])
        valid_pids = []
        for sib_id in valid_sib_ids:
            pdg_id = self.event_generator._lib.isib_pid2pdg(sib_id)
            valid_pids.append(pdg_id)
        
        return set(valid_pids) 
    
    def set_average_A(self, A):
        self.average_A = A
        self.factor_sigma = self.average_A * self.mass_barn

    def _sigma_wrapper(self):
        k = self.event_generator.event_kinematics

        if k.p1pdg not in self.valid_pids:
            raise Exception("Not a valid pdg for beam")

        sigma_hair = self.event_generator._lib.sib_sigma_hair

        kabs = abs(k.p1pdg)
        if 1000 < kabs < 10000:
            sigproj = 1
        elif kabs % 1000 < 300:
            sigproj = 2
        else:
            sigproj = 3
        sigma = sigma_hair(sigproj, k.ecm)
        if isinstance(sigma, tuple):
            return sigma[0]
        return sigma  
    
    
    def get_average_xdepth(self):
        self.event_generator.event_kinematics = impy.kinematics.FixedTarget(
            self.particle.energy, int(self.particle.pid), (14, 7)
        )
        return self.factor_sigma / self._sigma_wrapper()       
          

class DecayEvent(ParticleEvent):
    from particle_decay import MCEqDBFactory, ParticleDecay
    
    def __init__(self, particle):
        self.pdecay =  self.ParticleDecay(self.MCEqDBFactory().get_db(),
                                     particle.pid, particle.energy)
        
        self.xh_conversion = ConvertX2H()
        super().__init__(particle)
        
    
    def set_particle(self, particle):
        self.particle = particle
        self.pdecay.set_decayed_particle(self.particle.pid, self.particle.energy)
        
        
    def get_xdepth(self):
        
        
    @abstractmethod    
    def get_products(self):
        pass
    

    
                    
        
        
class CascadeEvent:
    mbarn_in_cm2 = 1e-27
    proton_mass_g = 1.672621e-24
    mass_barn = proton_mass_g / mbarn_in_cm2

    def __init__(self, emin_threshold, decay_particle):

        self.decay_particle = decay_particle
        self.xh_conversion = ConvertX2H()

        self.emin_threshold = emin_threshold

        ekin = impy.kinematics.FixedTarget(10, "proton", (14, 7))
        self.event_generator = impy.models.Sibyll23d(ekin)

        self.set_average_A(14)

        valid_sib_ids = [7, 8, 9, 10, 11, 12, 13, 14, -13, -14]
        valid_sib_ids.extend([34, 35, 36, 37, 38, 39])
        valid_sib_ids.extend([59, 60, 71, 72, 74, 75])
        valid_sib_ids.extend([87, 88, 89, 99, 27])
        valid_pids = []
        for sib_id in valid_sib_ids:
            pdg_id = self.event_generator._lib.isib_pid2pdg(sib_id)
            valid_pids.append(pdg_id)

        self.valid_pids = set(valid_pids)

    def set_average_A(self, A):
        self.average_A = A
        self.factor_sigma = self.average_A * self.mass_barn

    def _sigma_wrapper(self):
        k = self.event_generator.event_kinematics

        if k.p1pdg not in self.valid_pids:
            raise Exception("Not a valid pdg for beam")

        sigma_hair = self.event_generator._lib.sib_sigma_hair

        kabs = abs(k.p1pdg)
        if 1000 < kabs < 10000:
            sigproj = 1
        elif kabs % 1000 < 300:
            sigproj = 2
        else:
            sigproj = 3
        sigma = sigma_hair(sigproj, k.ecm)
        if isinstance(sigma, tuple):
            return sigma[0]
        return sigma

    def get_average_xlength(self, pid, energy):
        self.event_generator.event_kinematics = impy.kinematics.FixedTarget(
            energy, int(pid), (14, 7)
        )
        sigma = 1 / self._sigma_wrapper()
        return self.factor_sigma * sigma

    def get_xlength(self, pid, energy):
        xlen = self.get_average_xlength(pid, energy)
        return -np.log(1 - random.random()) * xlen

    def get_generation(self, shower_particle):
        if shower_particle.energy < self.emin_threshold:
            return None

        try:
            xlength = self.get_xlength(shower_particle.pid, shower_particle.energy)
        except Exception:
            return None

        if not xlength or math.isnan(xlength):
            return None

        xlength += shower_particle.xlength
        event = next(self.event_generator(1)).final_state()

        inter_particles = []
        final_particles = []

        for i in range(len(event)):
            ev = event[i]
            cas_prt = CascadeParticle(ev.pid, ev.en, xlength)
            if (ev.en < self.emin_threshold) or (ev.pid not in self.valid_pids):
                final_particles.append(cas_prt)
            else:
                inter_particles.append(cas_prt)

        return inter_particles, final_particles

    def get_interaction_length(self, shower_particle):
        try:
            xlength = self.get_xlength(shower_particle.pid, shower_particle.energy)
        except Exception:
            xlength = None

        if not xlength or math.isnan(xlength):
            xlength = None
            
        return xlength    
        
    def get_decay_length(self, shower_particle):
        self.decay_particle.set_decayed_particle(
            shower_particle.pid, shower_particle.energy
        )
        length = self.decay_particle.get_decay_length()
        
        if length:
            return self.xh_conversion.get_delta_xdepth(shower_particle.xlength, length)
        else:
            return None     
    
    
    def get_event_particles(self, shower_particle):
        
            
    
    def get_interaction_event(self, shower_particle, xlength):
        xlength += shower_particle.xlength
        event = next(self.event_generator(1)).final_state()

        inter_particles = []
        final_particles = []

        for i in range(len(event)):
            ev = event[i]
            cas_prt = CascadeParticle(ev.pid, ev.en, xlength)
            if (ev.en < self.emin_threshold) or (ev.pid not in self.valid_pids):
                final_particles.append(cas_prt)
            else:
                inter_particles.append(cas_prt)

        return inter_particles, final_particles

    def get_decay_event(self, shower_particle, xlength):

        self.decay_particle.set_decayed_particle(
            shower_particle.pid, shower_particle.energy
        )
        
        inter_particles = []
        final_particles = []
        
        for product in self.decay_particle.get_decay_products():
            p_pid = product[0]
            p_energy = product[1]
            cas_prt = CascadeParticle(p_pid, p_energy, xlength)
            if (p_energy < self.emin_threshold) or (p_pid not in self.valid_pids):
                final_particles.append(cas_prt)
            else:
                inter_particles.append(cas_prt)
        return inter_particles, final_particles            

    def interaction_or_decay_event(self, shower_particle, xlength):
        interaction_length = self.get_interaction_length(shower_particle, xlength)
        decay_length = self.get_decay_length(shower_particle, xlength)

        if decay_length:
            if interaction_length < decay_length:
                return self.get_interaction_event(shower_particle, interaction_length)
            else:
                return self.get_decay_event(shower_particle, decay_length)


class CascadeDriver:
    def __init__(self, shower_generator):
        self.shower_generator = shower_generator
        self.iterations = 0
        self.final_products = []

    def run(self, pid, energy):
        initial_particle = CascadeParticle(pid, energy, 0)
        pstack = [initial_particle]
        while pstack:
            self.iterations += 1
            # if self.iterations % 100000:
            #     print(f"Pstack = {len(pstack)}, Fstack = {len(self.final_products)},"
            #           f" Iterations = {self.iterations}\r")

            cur_particle = pstack.pop()
            current_generation = self.shower_generator.get_generation(cur_particle)
            if current_generation:
                pstack.extend(current_generation[0])
                self.final_products.extend(current_generation[1])
            else:
                self.final_products.append(cur_particle)

    def get_particles(self):
        return self.final_products

    def get_iterations(self):
        return self.iterations


class ConvertX2H:
    from MCEq.geometry.density_profiles import CorsikaAtmosphere

    length_units = {"cm": 1, "m": 1e2, "km": 1e5}

    def __init__(self):
        self.cka_obj = self.CorsikaAtmosphere("SouthPole", "December")
        self.cka_obj.set_theta(0.0)
        self.length_unit = self.length_units["cm"]

    def set_theta(self, theta):
        self.cka_obj.set_theta(theta)

    def set_length_unit(self, unit):
        self.length_unit = self.length_units[unit]

    def __call__(self, x):
        return self.cka_obj.X2h(x) / self.length_unit

    def get_length(self, x1, x2):
        """Return a length between atmospheric depth x1 and x2
        x1 < x2. Because X2h coverts only to height, zenith angle 
        is also taken into account.

        Args:
            x1 (float): _description_
            x2 (float): _description_

        Returns:
            float: length in a set length units
        """
        return (self.cka_obj.X2h(x2) - self.cka_obj.X2h(x1)) / (
            np.cos(self.cka_obj.thrad) * self.length_unit
        )
        
    def get_delta_xdepth(self, x1, length):
        """Returns delta x = x2 - x1 for the path starting at point
        x1 and run `length` cm

        Args:
            x1 (_type_): _description_
            length (_type_): _description_

        Returns:
            _type_: _description_
        """
        h1 = self.cka_obj.X2h(x1)
        h2 = h1 - length * np.cos(self.cka_obj.thrad)
        x2 = self.cka_obj.h2X(h2)
        return x2 - x1
