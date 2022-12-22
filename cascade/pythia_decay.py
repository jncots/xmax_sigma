import impy
from impy.common import EventData
from impy.util import mass, energy2momentum
import dataclasses
import numpy as np
from particle import Particle, PDGID


@dataclasses.dataclass
class SimpleEvent:
    pid = []
    status = []
    px = []
    py = []
    pz = []
    en = []
    m = []
    tau = []

    def append(self, pid, status, px, py, pz, en, m, tau):
        self.pid.append(pid)
        self.status.append(status)
        self.px.append(px)
        self.py.append(py)
        self.pz.append(pz)
        self.en.append(en)
        self.m.append(m)
        self.tau.append(tau)

    def append_1dz(self, pid, energy, *, status=1):
        self.pid.append(pid)
        self.status.append(status)
        self.en.append(energy)
        m = mass(pid)
        self.m.append(m)
        self.pz.append(energy2momentum(energy, m))
        self.px.append(0)
        self.py.append(0)
        self.tau.append(0)

    def __len__(self):
        return len(self.pid)

    def slice_it(self, arg):
        evt = SimpleEvent()
        evt.pid = self.pid[arg]
        evt.status = self.status[arg]
        evt.px = self.px[arg]
        evt.py = self.py[arg]
        evt.pz = self.pz[arg]
        evt.en = self.en[arg]
        evt.m = self.m[arg]
        evt.tau = self.tau[arg]
        return evt

    def pid_en_view(self):
        view = []
        for i in range(len(self.pid)):
            view.append((self.pid[i], self.en[i]))

        return view


class Pythia8DecayAfterburner:
    def __init__(self, stable_list):
        self.stable_list = stable_list
        self._init_pythia()

    def _init_pythia(self):

        ekin = impy.kinematics.CenterOfMass(1 * impy.constants.TeV, "proton", "proton")
        generator = impy.models.Pythia8(ekin)
        self.pythia = generator._lib.pythia
        self.pythia.readString("ProcessLevel:all = off")
        self.pythia.readString("ParticleDecays:tau0Max = 1e100")
        # Set muons unstable
        # self.pythia.particleData.mayDecay(13, True)
        self.pythia.init()
        print(f"Finish of initialization")

    def __call__(self, event):

        init_size = len(event)

        for ip in range(len(event)):
            if event.status[ip] != 1 or abs(event.pid[ip]) in self.stable_list:
                continue

            # Set particle to not final and simulate decay
            event.status[ip] = 2
            self.pythia.event.reset()
            # put decaying particle
            self.pythia.event.append(
                int(event.pid[ip]),
                91,
                0,
                0,
                event.px[ip],
                event.py[ip],
                event.pz[ip],
                event.en[ip],
                event.m[ip],
            )
            self.pythia.particleData.mayDecay(int(event.pid[ip]), True)
            # Decay it
            self.pythia.forceHadronLevel()
            for p in self.pythia.event:
                if not p.isFinal():
                    continue
                event.append(p.id(), 1, p.px(), p.py(), p.pz(), p.e(), p.m(), p.tau())

        return event.slice_it(slice(init_size, None, None))


class DecayByPythia:
    def __init__(self):
        self.afterburner = Pythia8DecayAfterburner(
            stable_list=[2212, 11, 12, 14, 15, 16, 22]
        )

    def get_decayed_products(self, pid, energy):
        event = SimpleEvent()
        event.append_1dz(pid, energy)
        event = self.afterburner(event)
        return event.pid_en_view()



# decay_event = DecayByPythia()

# print(decay_event.get_decayed_products(211, 1e3))

# event = SimpleEvent()
# event.append_1dz(3222, 1e3)
# # event.append_1dz(-211, 1e3)
# # event.append_1dz(2112, 1e3)
# # event.append_1dz(11, 1e1)
# # event.append_1dz(22, 1e2)

# # event = event[1:]

# afterburner = Pythia8DecayAfterburner(stable_list=[2212, 11, 12, 14, 15, 16, 22])
# event = afterburner(event)

# for ip in range(len(event)):
#     print(
#         f"id = {ip}, pid = {event.pid[ip]}, "
#         f"status = {event.status[ip]}, "
#         f"energy = {event.en[ip]}, "
#         f"mass = {event.m[ip]}, "
#         f"tau = {event.tau[ip]}"
#     )

# print(event.pid_en_view())


# pythia.readString("ProcessLevel:all = off")
# pythia.readString("ParticleDecays:tau0Max = 1e100")
#         # Set muons unstable
# # pythia.particleData.mayDecay(13, True)
# # pythia.init()

# # pythia.event.reset()

# event = SimpleEvent()
# event.append_1dz(211, 1e3)
# ip = 0

# pythia.event.append(
#     int(event.pid[ip]),
#     91,
#     0,
#     0,
#     event.px[ip],
#     event.py[ip],
#     event.pz[ip],
#     event.en[ip],
#     event.m[ip],
# )

# pythia.particleData.mayDecay(int(event.pid[ip]), True)
# # Decay it
# pythia.forceHadronLevel()
# for p in pythia.event:
#     # if not p.isFinal():
#     #     continue
#     event.append(p.id(), p.status(), p.px(), p.py(), p.pz(), p.e(), p.m(), p.tau())


# for ip in range(len(event)):
#     print(f"id = {ip}, pid = {event.pid[ip]}, "
#           f"status = {event.status[ip]}, "
#           f"energy = {event.en[ip]}, "
#           f"mass = {event.m[ip]}, "
#           f"tau = {event.tau[ip]}")


# print(pythia.particleData.findParticle(13).m0)
# # for event in generator(10):
# #     print(event.pid)
