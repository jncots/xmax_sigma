from particle_event import ParticleEvent, CascadeParticle
from xdepth_conversion import XdepthConversion
from particle_decay import MCEqDBFactory, ParticleDecay


class DecayEvent(ParticleEvent):
    def __init__(self, particle):
        self.pdecay = ParticleDecay(
            MCEqDBFactory().get_db(), particle.pid, particle.energy
        )

        self._get_prod_ncalls = 0
        self.xd_convertion = XdepthConversion()
        super().__init__(particle)

    def set_particle(self, particle):
        self.particle = particle
        self.pdecay.set_decayed_particle(self.particle.pid, self.particle.energy)

    def get_xdepth(self):
        self.last_xdepth = None

        length = self.pdecay.get_decay_length()
        if not length or length > 200e5:
            return None


        self.last_xdepth = self.xd_convertion.get_delta_xdepth(
            self.particle.xdepth, length
        )
   

        return self.last_xdepth

    def get_products(self):

        if not self.last_xdepth:
            return None

        total_xdepth = self.particle.xdepth + self.last_xdepth

        if total_xdepth > self.max_xdepth:
            return None

        event = self.pdecay.get_decay_products()

        generation_number = self.particle.generation_number + 1
        products = []
        for particle in event:
            cp = CascadeParticle(
                particle[0], particle[1], total_xdepth, 2, generation_number
            )
            cp.parent.append(self.particle)
            products.append(cp)

        self.last_xdepth = None

        self._get_prod_ncalls += 1
        return products



de = DecayEvent(CascadeParticle(-13, 100, 100))

print(de.pdecay.get_decay_length())