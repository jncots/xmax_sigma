import impy
import random
import numpy as np
import dataclasses
import math


@dataclasses.dataclass
class ShowerParticle:
    id: int
    parent_id: int
    interacted: bool
    pid: int
    energy: float
    xlength: float


class ShowerGenerator:
    emin_threshold = 1e4
    mbarn_in_cm2 = 1e-27
    proton_mass_g = 1.672621e-24
    mass_barn = proton_mass_g / mbarn_in_cm2
    id_in_shower = 0

    event_generator = None

    def __init__(self, emin_threshold):
        print(f"Evgen = {ShowerGenerator.event_generator}")
        if self.event_generator is None:
            ekin = impy.kinematics.FixedTarget(10, "proton", (14, 7))
            ShowerGenerator.event_generator = impy.models.Sibyll23d(ekin)
            self._set_average_A(14)
            self.event_generator_is_set = True
        valid_sib_ids = [7, 8, 9, 10, 11, 12, 13, 14, -13, -14]
        valid_sib_ids.extend([34, 35, 36, 37, 38, 39])
        valid_sib_ids.extend([59, 60, 71, 72, 74, 75])
        valid_sib_ids.extend([87, 88, 89, 99, 27])
        self.valid_sib_ids = set(valid_sib_ids)

    def _set_average_A(self, A):
        self.average_A = A
        self.factor_sigma = A * self.mass_barn

    def _sigma_wrapper(self):
        k = ShowerGenerator.event_generator.event_kinematics

        sib_id = ShowerGenerator.event_generator._lib.isib_pdg2pid(k.p1pdg)
        if sib_id not in self.valid_sib_ids:
            raise Exception("Not a valid sib_id")

        sigma_hair = ShowerGenerator.event_generator._lib.sib_sigma_hair

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
        ShowerGenerator.event_generator.event_kinematics = impy.kinematics.FixedTarget(
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

        generation = []

        for i in range(len(event)):
            ev = event[i]
            # The id_in_shower for initial particle is 0
            # So the last id_in_shower is also number of generated particles so far
            self.id_in_shower += 1
            sp = ShowerParticle(
                self.id_in_shower, shower_particle.id, False, ev.pid, ev.en, xlength
            )
            generation.append(sp)

        return generation


class CascadeEngine:
    def __init__(self, shower_generator):
        self.shower_generator = shower_generator
        self.iterations = 0
        self.interactions = 0
        self.final_products = []

    def run(self, pid, energy):
        initial_particle = ShowerParticle(0, 0, False, pid, energy, 0)
        pstack = [initial_particle]
        while pstack:
            self.iterations += 1
            # if self.iterations % 10000:
            #     print(f"Pstack = {len(pstack)}, Fstack = {len(self.final_products)},"
            #           f" Iterations = {self.iterations}\r")
            
            cur_particle = pstack.pop()
            current_generation = self.shower_generator.get_generation(cur_particle)
            if current_generation:
                self.interactions += 1
                pstack.extend(current_generation)
            else:
                self.final_products.append(cur_particle)

        return self.final_products

    def get_shower(self):
        return self.final_products

    def get_iterations(self):
        return self.iterations

    def get_interactions(self):
        return self.interactions


class ConvertX2H:
    from MCEq.geometry.density_profiles import CorsikaAtmosphere

    def __init__(self):
        self.cka_obj = self.CorsikaAtmosphere("SouthPole", "December")
        self.cka_obj.set_theta(0.0)

    def set_theta(self, theta):
        self.cka_obj.set_theta(theta)

    def __call__(self, x):
        return self.cka_obj.X2h(x) / 1e5
