import impy
import random
import numpy as np
import dataclasses
import math


@dataclasses.dataclass
class CascadeParticle:
    id: int
    parent_id: int
    interacted: bool
    pid: int
    energy: float
    xlength: float


class CascadeEvent:
    mbarn_in_cm2 = 1e-27
    proton_mass_g = 1.672621e-24
    mass_barn = proton_mass_g / mbarn_in_cm2

    def __init__(self, emin_threshold):
        
        self.emin_threshold = emin_threshold
        self.id_in_shower = 0
        
        ekin = impy.kinematics.FixedTarget(10, "proton", (14, 7))
        self.event_generator = impy.models.Sibyll23d(ekin)
        
        self._set_average_A(14)
            
        valid_sib_ids = [7, 8, 9, 10, 11, 12, 13, 14, -13, -14]
        valid_sib_ids.extend([34, 35, 36, 37, 38, 39])
        valid_sib_ids.extend([59, 60, 71, 72, 74, 75])
        valid_sib_ids.extend([87, 88, 89, 99, 27])
        valid_pids = []
        for sib_id in valid_sib_ids:
            pdg_id = self.event_generator._lib.isib_pid2pdg(sib_id)
            valid_pids.append(pdg_id)
        
        self.valid_pids = set(valid_pids)

    def _set_average_A(self, A):
        self.average_A = A
        self.factor_sigma = A * self.mass_barn

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

        # Means that we used the particle
        # Do we need this ?
        shower_particle.interacted = True

        if shower_particle.energy < self.emin_threshold:
            return None

        try:
            xlength = self.get_xlength(shower_particle.pid, shower_particle.energy)
        except Exception:
            return None

        if not xlength or math.isnan(xlength):
            return None

        xlength += shower_particle.xlength
        # print(f"Shower pid {shower_particle.pid}")
        event = next(self.event_generator(1)).final_state()

        inter_particles = []
        final_particles = []
        

        for i in range(len(event)):
            ev = event[i]
            # The id_in_shower for initial particle is 0
            # So the last id_in_shower is also number of generated particles so far
            self.id_in_shower += 1
            cas_prt = CascadeParticle(
                self.id_in_shower, shower_particle.id, False, ev.pid, ev.en, xlength
            )
            if (ev.en < self.emin_threshold) or (ev.pid not in self.valid_pids):
                final_particles.append(cas_prt)
            else:
                inter_particles.append(cas_prt)   

        return inter_particles, final_particles


class CascadeDriver:
    def __init__(self, shower_generator):
        self.shower_generator = shower_generator
        self.iterations = 0
        self.final_products = []

    def run(self, pid, energy):
        initial_particle = CascadeParticle(0, 0, False, pid, energy, 0)
        pstack = [initial_particle]
        while pstack:
            self.iterations += 1
            # if self.iterations % 10000:
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

    def __init__(self):
        self.cka_obj = self.CorsikaAtmosphere("SouthPole", "December")
        self.cka_obj.set_theta(0.0)

    def set_theta(self, theta):
        self.cka_obj.set_theta(theta)

    def __call__(self, x):
        return self.cka_obj.X2h(x) / 1e5
