from particle import Particle
import numpy as np
import math
from particletools.tables import PYTHIAParticleData


class TabParticles:
    """Tabulate properties of particles as numpy arrays
    for fast access. The following properties are tabulated

    TabParticles.mass(pdg_array)
    TabParticles.ctau(pdg_array)
    """

    # Absolute pdg id which is tabulated in numpy array
    # Properties of particles with abs(pdg) > _pdg_max
    # are not be tabulated and are obtained from
    # original function, i.e. they work slower
    _pdg_max = 5000
    _none_value = np.int16(32_767)
    _GeV = 1e-3
    _cm = 1e-1

    def __init__(self):
        self._init_all_particles()
        self.pythia_pdata = PYTHIAParticleData()
        self._construct_maps()
        self._fill_mass_tab()
        self._fill_ctau_tab()

    def _init_all_particles(self):
        self._all_particles_list = Particle.findall()
        self._all_particles_dict = {p.pdgid: p for p in self._all_particles_list}

    def _construct_maps(self):
        self._pdg2id = np.empty(2 * self._pdg_max + 1, dtype=np.int16)
        self._pdg2id.fill(self._none_value)

        id2pdg = np.empty(2 * self._pdg_max + 1, dtype=np.int32)

        pid = 0
        for p in self._all_particles_list:
            pdg = int(p.pdgid)
            if abs(pdg) <= self._pdg_max:
                self._pdg2id[pdg] = pid
                id2pdg[pid] = pdg
                pid += 1

        self._id2pdg = np.empty(pid, dtype=np.int32)
        self._id2pdg[:] = id2pdg[0:pid]
        self._len_tab = len(self._id2pdg)

    def _get_ctau(self, pdg):
        """Returns c*tau, in cm"""
        ctau = self._all_particles_dict[pdg].ctau
        if (ctau is None) or math.isinf(ctau):
            return np.float64(0)
        else:
            return np.float64(ctau * self._cm)

    def _fill_ctau_tab(self):
        self._ctau = np.zeros(self._len_tab, dtype=np.float64)
        for pid, pdg in enumerate(self._id2pdg):
            self._ctau[pid] = self._get_ctau(pdg)

    def _get_mass(self, pdg):
        """Returns mass in GeV"""
        if self._all_particles_dict[pdg].mass is None:
            return np.float64(0)
        else:
            return np.float64(self._all_particles_dict[pdg].mass) * self._GeV

    def _fill_mass_tab(self):
        self._mass = np.zeros(self._len_tab, dtype=np.float64)
        for pid, pdg in enumerate(self._id2pdg):
            self._mass[pid] = self._get_mass(pdg)

    def mass(self, pdg_array):
        try:
            mass_array = self._mass[self._pdg2id[pdg_array]]
        except (IndexError):
            # It is assumed to be a rare case
            mass_array = np.frompyfunc(self._get_mass, 1, 1)(pdg_array).astype(
                "float64"
            )
        return mass_array

    def ctau(self, pdg_array):
        try:
            ctau_array = self._ctau[self._pdg2id[pdg_array]]
        except (IndexError):
            # It is assumed to be a rare case
            ctau_array = np.frompyfunc(self._get_ctau, 1, 1)(pdg_array).astype(
                "float64"
            )
        return ctau_array


if __name__ == "__main__":
    tp = TabParticles()

    print(tp.ctau([111, 11, 13, -13]))
    print(tp.ctau([-1001112720, 1000902270, 1000721650]))

    print(tp.mass(np.array([111, 111, 111])))
    print(tp.mass([-1001112720, 1000902270, 1000721650]))
