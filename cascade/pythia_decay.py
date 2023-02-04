import chromo
import numpy as np
from pathlib import Path
from particle_event import CascadeParticle

chormo_path = Path(chromo.__file__).parent


class Pythia8DecayAfterburner:
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

    def _newcall(self, event):

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

    def get_number_of_decays(self):
        return self.number_of_decays

from contextlib import contextmanager
@contextmanager
def _send_to_null():
    """Send stdout output to devnull
    
    Usage:
        with _send_to_null():
            ...code to quite down
    """
    # Solution from https://stackoverflow.com/a/6735958
    import os
    import sys
        
    original_stdout = sys.stdout
    fnull = open(os.devnull, 'w')

    sys.stdout = fnull  
    yield
    sys.stdout = original_stdout   

import chromo
class Pythia8Decay(chromo.models.Pythia8):
    def __init__(self, *, seed=None):
        from random import randint

        if seed is None:
            seed = randint(1, 10000000)

        evt_kin = chromo.kinematics.CenterOfMass(1000, 2212, 2212)
        with _send_to_null():
            super().__init__(evt_kin, seed=seed)

    def _set_kinematics(self, kin):
        """Overriding _set_kinematics for decays only
        """
        pythia = self._lib.pythia
        pythia.settings.resetAll()

        config = [
            "Random:setSeed = on",
            f"Random:seed = {self._seed}",
            "Print:quiet = on",
            "ProcessLevel:all = off",
            "ParticleDecays:tau0Max = 1e100",
        ]

        for line in config:
            if not pythia.readString(line):
                raise RuntimeError(f"readString({line!r}) failed")

        if not pythia.init():
            raise RuntimeError("Pythia8 initialization failed")
        
        self.pythia = pythia
        
    def run_decay(self, event):
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


if __name__ == "__main__":

    # ff = Pythia8Decay()
       

    # afterburner = Pythia8DecayAfterburner()
    # event111 = [
    #     CascadeParticle(2212, 1e3, 77),
    #     CascadeParticle(111, 13e2, 888)
    # ]

    # event = afterburner._newcall(event111)
    # # event = afterburner(event111)

    # for p in event:
    #     print(p)
