import chromo
from chromo.util import energy2momentum
from pathlib import Path
from particle_event import CascadeParticle

chormo_path = Path(chromo.__file__).parent

class Pythia8DecayAfterburner:
    def __init__(self):
        self._init_pythia()

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
        
        new_event = []
        
        self.pythia.event.reset()
        
        for particle in event:
            mp = self.pythia.particleData.findParticle(particle.pid).m0
            self.pythia.event.append(
                particle.pid,
                91,
                0,
                0,
                0,
                0,
                energy2momentum(particle.energy, mp),
                particle.energy,
                mp,
            )
            self.pythia.particleData.mayDecay(particle.pid, True)
            
        # Decay it
        self.pythia.forceHadronLevel()
        ievent = 0
        for p in self.pythia.event:
            if ievent == 0:
                ievent += 1
                continue
        
            if p.mother1() == 0:
                if p.isFinal():
                    new_event.append(event[ievent - 1])
            else:
                mother = event[p.mother1() - 1]
                cp = CascadeParticle(pid = p.id(), 
                                    energy = p.e(), 
                                    xdepth = mother.xdepth, 
                                    production_mode = 2,
                                    generation_number = mother.generation_number + 1,
                                    parent = [mother])    
                new_event.append(cp)
            ievent += 1
        return new_event


