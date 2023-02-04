import random
import numpy as np
import math

import chromo
from particle_event import CascadeParticle
from xdepth_conversion import XdepthConversion


class HadronEvent:
    mbarn_in_cm2 = 1e-27
    proton_mass_g = 1.672621e-24
    mass_barn = proton_mass_g / mbarn_in_cm2

    def __init__(self):
        ekin = chromo.kinematics.FixedTarget(20000, "proton", (14, 7))
        self.event_generator = chromo.models.Sibyll23d(ekin)

        self.set_average_A(14)
        self.set_target((14, 7))
        self.valid_pids = self._get_valid_pids()

        self.xdepth_max = XdepthConversion().get_max_xdepth()
        self._get_prod_ncalls = 0
        self.valid_pids = self._get_valid_pids()

    def set_average_A(self, A):
        self.average_A = A
        self.factor_sigma = self.average_A * self.mass_barn

    def set_target(self, target):
        self.target = target

    def get_xdepth(self, particle):

        if particle.pid not in self.valid_pids:
            return None

        try:
            average_xdepth = self.get_average_xdepth(particle)
        except Exception:
            return None

        if not average_xdepth or math.isnan(average_xdepth):
            return None

        return -np.log(1 - random.random()) * average_xdepth

    def get_products(self, particle, xdepth):        
        
        total_xdepth = particle.xdepth + xdepth
        if total_xdepth > self.xdepth_max:
            particle.xdepth = self.xdepth_max
            particle.final_code = 3
            return None

        self.event_generator.kinematics = chromo.kinematics.FixedTarget(
            particle.energy, particle.pid, self.target
        )
        event = next(self.event_generator(1)).final_state()

        generation_number = particle.generation_number + 1
        products = []
        for i in range(len(event.pid)):
            cp = CascadeParticle(
                int(event.pid[i]), event.en[i], total_xdepth, 1, generation_number
            )
            cp.parent.append(particle)
            products.append(cp)
        self._get_prod_ncalls += 1
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

    def _sigma_wrapper(self, pid, energy):
        event_kin = chromo.kinematics.FixedTarget(energy, pid, self.target)

        if int(event_kin.p1) not in self.valid_pids:
            print(f"{event_kin.p1} is not in valid_pids")
            raise Exception("Not a valid pdg for beam")

        self.event_generator.kinematics = event_kin
        sigma_hair = self.event_generator._lib.sib_sigma_hair

        pid_abs = abs(int(event_kin.p1))
        if 1000 < pid_abs < 10000:
            sigproj = 1
        elif pid_abs % 1000 < 300:
            sigproj = 2
        else:
            sigproj = 3
        sigma = sigma_hair(sigproj, event_kin.ecm)
        if isinstance(sigma, tuple):
            return sigma[0]
        return sigma

    def get_average_xdepth(self, particle):
        return self.factor_sigma / self._sigma_wrapper(particle.pid, particle.energy)


if __name__ == "__main__":
    hev = HadronEvent()
    particle = CascadeParticle(2212, 30000, 0.0)

    print(f"xdepth_average = {hev.get_average_xdepth(particle)}")
    for i in range(10):
        print(f"i = {i}, xdepth = {hev.get_xdepth(particle)}")

    xdepth = hev.get_xdepth(particle)
    products = hev.get_products(particle, xdepth)
    print(products)
