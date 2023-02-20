import chromo
import numpy as np
from pathlib import Path
# from particle_event import CascadeParticle

chormo_path = Path(chromo.__file__).parent


# class Pythia8DecayAfterburner:
#     def __init__(self):
#         self._init_pythia()
#         self.number_of_decays = 0

#     def _init_pythia(self):
#         import importlib
#         from random import randint

#         lib = importlib.import_module(f"chromo.models._pythia8")
#         xml_path = chormo_path / "iamdata/Pythia8/xmldoc"
#         self.pythia = lib.Pythia(str(xml_path), False)
#         seed = randint(1, 10000000)
#         self.pythia.settings.resetAll()
#         self.pythia.readString("Random:setSeed = on")
#         self.pythia.readString(f"Random:seed = {seed}")
#         self.pythia.readString("Print:quiet = on")
#         self.pythia.readString("ProcessLevel:all = off")
#         self.pythia.readString("ParticleDecays:tau0Max = 1e100")
#         self.pythia.init()




class DecayInteraction:
    def __init__(self):
        self._init_pythia()
        self.number_of_decays = 0
        
    def _init_pythia(self):
        import importlib
        from random import randint

        lib = importlib.import_module(f"chromo.models._pythia8")
        xml_path = chormo_path / "iamdata/Pythia8/xmldoc"
        self.pythia = lib.Pythia(str(xml_path), False)
        seed = randint(1, 10000000)
        self.pythia.settings.resetAll()
        self.pythia.readString("Random:setSeed = on")
        self.pythia.readString(f"Random:seed = {seed}")
        self.pythia.readString("Print:quiet = on")
        self.pythia.readString("ProcessLevel:all = off")
        self.pythia.readString("ParticleDecays:tau0Max = 1e100")
        self.pythia.init()
        
    def __call__(self, event):

        self.number_of_decays = 0
        all_particles = []
        final_particles = []

        self.pythia.event.reset()

        for particle in event:
            m0 = self.pythia.particleData.findParticle(particle.pid).m0
            self.pythia.event.append(
                particle.pid,
                91,
                0,
                0,
                0,
                0,
                np.sqrt((particle.energy + m0) * (particle.energy - m0)),  # pz
                particle.energy,
                m0,
            )
            self.pythia.particleData.mayDecay(particle.pid, True)

        # Decay it
        self.pythia.forceHadronLevel()
        ievent = 0
        for p in self.pythia.event:
            # The first record in self.pythia.event is
            # the event itself, not a particle, therefore pass
            if ievent == 0:
                # Add None for one to one mapping
                # So mother1() will be an index in all_particles
                all_particles.append(None)
                ievent += 1
                continue

            if ievent <= len(event):
                all_particles.append(event[ievent - 1])
            else:
                mother = all_particles[p.mother1()]
                particle = CascadeParticle(
                    pid=p.id(),
                    energy=p.e(),
                    xdepth=mother.xdepth_decay,
                    xdepth_decay=mother.xdepth_decay,
                    production_mode=2,
                    generation_number=mother.generation_number + 1,
                    parent=[mother],
                )
                all_particles.append(particle)

            if p.isFinal():
                if p.mother1() == 0:
                    final_particles.append(all_particles[ievent])
                else:
                    mother = all_particles[p.mother1()]
                    particle = CascadeParticle(
                        pid=p.id(),
                        energy=p.e(),
                        xdepth=mother.xdepth_decay,
                        xdepth_decay=mother.xdepth_decay,
                        production_mode=2,
                        generation_number=mother.generation_number + 1,
                        parent=[mother],
                    )
                    final_particles.append(particle)
            else:
                self.number_of_decays += 1

            ievent += 1
        return final_particles        
    
    

if __name__ == "__main__":
    
    dint = DecayInteraction()    