from particle_array import ParticleArray
from tab_pproperties import TabulatedParticleProperties
from xdepth_on_table import XdepthOnTable
import numpy as np


class DecayXdepth:
    def __init__(self, 
                 particle_propeties = TabulatedParticleProperties(),
                 xdepth_on_table = XdepthOnTable()):
        self.pp_tab = particle_propeties
        self.xd_tab = xdepth_on_table

    def get_xdepth(self, pdg, energy, xdepth):
        mass = self.pp_tab.mass(pdg)
        gamma = np.divide(
            energy, mass, out=np.full_like(mass, np.inf), where=mass != 0
        )
        bgamma = np.sqrt((gamma + 1) * (gamma - 1))
        rnd = -np.log(1 - np.random.rand(len(pdg)))
        length = rnd * bgamma * self.pp_tab.ctau(pdg)
        print("length = ", length/1e5)
        xdepth[:] = self.xd_tab.add_len2x(xdepth, length)


if __name__ == "__main__":
    particle_stack = ParticleArray(size=100)
    decay_xdepth = DecayXdepth()

    particle_stack.push(
        pid=np.array([111, 22, 111, 13, -13, 2212]),
        energy=np.array([2e10, 2e0, 2e10, 1e0, 1e0, 1e0]),
        xdepth=np.array([400, 0, 0, 0, 0, 0]),
    )

    print(particle_stack.xdepth_decay[0 : len(particle_stack)])

    decay_xdepth.get_xdepth(particle_stack.pid[0 : len(particle_stack)],
        particle_stack.energy[0 : len(particle_stack)],
        particle_stack.xdepth_decay[0 : len(particle_stack)])

    print(particle_stack.xdepth_decay[0 : len(particle_stack)])