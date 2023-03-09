import chromo

evt_kin = chromo.kinematics.CenterOfMass(100, 2212, 2212)
gen = chromo.models.DpmjetIII191(evt_kin)

for event in gen(1):
    # self._lib.dtiont.lout
    # elf._lib.dtflka.lout
    print(f"FOUT = {gen._lib.dtflka.lout}")
    gen._lib.dt_evtout(4)
    print(event.pid)
    





# import numpy as np

# v1 = np.array([2, 5])
# v2 = np.array([3, 4])
# v3 = np.array([4, 3])
# v4 = np.array([5, 1])

# aa = np.array([v1, v2, v3, v4])
# bb = np.amin(aa, axis = 0)

# print(aa.min(0))
# print(aa.argmin(0))
