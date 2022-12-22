from cascade_event import CascadeEvent
from particle_event import CascadeParticle




class CascadeDriver:
    def __init__(self, cascade_event):
        self.cascade_event = cascade_event
        self.iterations = 0
        self.final_products = []

    def run(self, initial_particle):
        pstack = [initial_particle]
        while pstack:
            self.iterations += 1
            # if self.iterations % 100000:
            #     print(f"Pstack = {len(pstack)}, Fstack = {len(self.final_products)},"
            #           f" Iterations = {self.iterations}\r")

            cur_particle = pstack.pop()
            current_generation = self.cascade_event.get_event_particles(cur_particle)
            if current_generation:
                pstack.extend(current_generation[0])
                self.final_products.extend(current_generation[1])
            else:
                self.final_products.append(cur_particle)

    def get_particles(self):
        return self.final_products

    def get_iterations(self):
        return self.iterations