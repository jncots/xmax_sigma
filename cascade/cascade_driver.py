from cascade_event import CascadeEvent
from particle_event import CascadeParticle


class CascadeDriver:
    def __init__(self, cascade_event):
        self.cascade_event = cascade_event
        self.cascade_event._reset_ncalls()
        self.iterations = 0
        self.final_products = []

    def run(self, initial_particle):

        decaying_particles = []
        pstack = [initial_particle]
        while pstack:
            self.iterations += 1
            if self.iterations % 1000 == 0:
            # if True:
                print(
                    f"Pstack = {len(pstack)}, Fstack = {len(self.final_products)},"
                    f" Iterations = {self.iterations},"
                    f" Decaying = {len(decaying_particles)}\r"
                )

            cur_particle = pstack.pop()
            current_generation = self.cascade_event.get_event_particles(cur_particle)
            # print(f"Current_gen[0] = {len(current_generation[0])}")
            # print(f"Particles on stack = {len(pstack)}")
            if len(current_generation[0]) > 0:
                pstack.extend(current_generation[0])
            if len(current_generation[1]) > 0:
                self.final_products.extend(current_generation[1])
            if len(current_generation[2]) > 0:
                decaying_particles.extend(current_generation[2])

            if not pstack:
                # print(f"Pstack is empty")
                # input()
                if decaying_particles:
                    pstack.extend(self.cascade_event.afterburner(decaying_particles))
                    decaying_particles = []

        self._set_statistics_variables()

    def _set_statistics_variables(self):
        self.num_decays = self.cascade_event.decay_event._get_prod_ncalls
        self.num_interactions = self.cascade_event.hadron_event._get_prod_ncalls
        self.num_events = self.iterations
        self.num_final_particles = len(self.final_products)

    def get_particles(self):
        return self.final_products

    def get_iterations(self):
        return self.iterations


if __name__ == "__main__":
    cas_driver = CascadeDriver(CascadeEvent(1e3))
    particle = CascadeParticle(2212, 1e8, 0)

    cas_driver.run(particle)
