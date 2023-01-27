import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parents[1]))

from cascade.pythia_decay import Pythia8DecayAfterburner
from cascade.particle_event import CascadeParticle


afterburner = Pythia8DecayAfterburner()
event111 = [
    CascadeParticle(2212, 1e3, 77),
    CascadeParticle(111, 13e2, 888)
]

event = afterburner(event111)

for p in event:
    print(p)

