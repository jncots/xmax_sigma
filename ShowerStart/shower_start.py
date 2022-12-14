import impy
import random
import numpy as np
import dataclasses
import math


@dataclasses.dataclass
class CascadeParticle:
    pid: int
    energy: float
    xlength: float


class CascadeEvent:
    mbarn_in_cm2 = 1e-27
    proton_mass_g = 1.672621e-24
    mass_barn = proton_mass_g / mbarn_in_cm2

    def __init__(self, emin_threshold):

        self.emin_threshold = emin_threshold

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

    def interaction_or_decay_event():
        interaction_length = get_interaction_length()
        decay_length = get_decay_length()

        if decay_length:
            if interaction_length < decay_length:
                return interaction_length, get_interaction_event()
            else:
                return decay_length, get_decay_event()


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


class DecayLength:
    from particletools.tables import PYTHIAParticleData

    def __init__(self, pdg, energy):
        pythia_pdata = self.PYTHIAParticleData()
        self.c_decay_time = pythia_pdata.ctau(pdg)
        self.mass = pythia_pdata.mass(pdg)
        self.set_energy(energy)

    def set_energy(self, energy):
        gamma = energy / self.mass
        beta_gamma = np.sqrt((gamma + 1) * (gamma - 1))
        self.decay_length = beta_gamma * self.c_decay_time

    def _random_value(self):
        return -np.log(1 - random.random())

    def get_decay_length(self):
        return self.decay_length * self._random_value()


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
        return (self.cka_obj.X2h(x2) - self.cka_obj.X2h(x1)) / (
            np.cos(self.cka_obj.thrad) * self.length_unit
        )
