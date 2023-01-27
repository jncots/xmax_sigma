import random
import numpy as np
import math

import chromo
from chromo.kinematics import EventFrame
from particle_event import ParticleEvent, CascadeParticle


class HadronEvent(ParticleEvent):

    mbarn_in_cm2 = 1e-27
    proton_mass_g = 1.672621e-24
    mass_barn = proton_mass_g / mbarn_in_cm2

    def __init__(self, particle):
        super().__init__(particle)

        ekin = chromo.kinematics.FixedTarget(20000, "proton", (14, 7))
        self.event_generator = chromo.models.Sibyll23d(ekin)

        self._get_prod_ncalls = 0
        self.set_average_A(14)
        self.valid_pids = self._get_valid_pids()

    def set_particle(self, particle):
        self.particle = particle

    def get_xdepth(self):
        """Returns slant depth of next interaction or None"""

        self.last_xdepth = None

        if self.particle.pid not in self.valid_pids:
            return None

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

        if total_xdepth > self.max_xdepth:
            return None

        event = next(self.event_generator(1)).final_state()

        generation_number = self.particle.generation_number + 1
        products = []
        for i in range(len(event.pid)):
            cp = CascadeParticle(
                int(event.pid[i]), event.en[i], total_xdepth, 1, generation_number
            )
            cp.parent.append(self.particle)
            products.append(cp)

        self.last_xdepth = None
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

    def set_average_A(self, A):
        self.average_A = A
        self.factor_sigma = self.average_A * self.mass_barn

    def _sigma_wrapper(self):
        k = self.event_generator.kinematics

        if int(k.p1) not in self.valid_pids:
            raise Exception("Not a valid pdg for beam")

        sigma_hair = self.event_generator._lib.sib_sigma_hair

        kabs = abs(int(k.p1))
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
        self.event_generator.kinematics = chromo.kinematics.FixedTarget(
            self.particle.energy, int(self.particle.pid), (14, 7)
        )
        return self.factor_sigma / self._sigma_wrapper()
