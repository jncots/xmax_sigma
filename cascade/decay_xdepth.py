from particle_array import ParticleArray
from tab_particles import TabParticles
from tab_xdepth import TabXdepth
import numpy as np


class DecayXdepth:
    def __init__(self, theta_deg=0):
        self.tab_pt = TabParticles()
        self.tab_xd = TabXdepth(theta_deg)

    def set_xdepth_decay(self, pstack):

        pslice = slice(0, len(pstack))
        pstack.energy[pslice]

        mass = self.tab_pt.mass(pstack.pid[pslice])
        gamma = np.divide(
            pstack.energy[pslice], mass, out=np.full_like(mass, np.inf), where=mass != 0
        )

        bgamma = np.sqrt((gamma + 1) * (gamma - 1))
        rnd = -np.log(1 - np.random.rand(len(pstack)))
        length = rnd * bgamma * self.tab_pt.ctau(pstack.pid[pslice])
        # print(f"ctau = {self.tab_pt.ctau(pstack.pid[pslice])}")
        # print(f"length = {length}")
        # print(f"bgamma = {bgamma}")
        pstack.xdepth_decay[pslice] = self.tab_xd.add_len2x(
            pstack.xdepth[pslice], length
        )


if __name__ == "__main__":
    particle_stack = ParticleArray(size=100)
    decay_xdepth = DecayXdepth(theta_deg=0)

    particle_stack.push(
        pid=np.array([111, 22, 111, 13, -13, 2212]),
        energy=np.array([2e10, 2e3, 2e10, 1e0, 1e0, 1e0]),
        xdepth=np.array([0, 0, 0, 0, 0, 0]),
    )

    print(particle_stack.xdepth_decay[0 : len(particle_stack)])

    decay_xdepth.set_xdepth_decay(particle_stack)

    print(particle_stack.xdepth_decay[0 : len(particle_stack)])
