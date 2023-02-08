import chromo
from particle import PDGID, Particle
import numpy as np


all_particles = Particle.findall()
# # print(len(all_particles))

pd = dict()
ppp = dict()
pnum = dict()
ii = 0

# print(all_particles)

# for p in all_particles:  
# print([(int(p.pdgid), p.name) for p in all_particles])

    
    
for p in all_particles:
    # hhash = int(p.pdgid) % 6521711
    
    # print(int(p.pdgid))
    # # print(int(p.pdgid), np.binary_repr(int(p.pdgid)))
    # input()

    if int(p.pdgid) > 0 and int(p.pdgid) < 30000:
        hhash = int(p.pdgid) % 5273

        pd[hhash] = pd.get(hhash, 0) + 1
        ppp[hhash] = int(p.pdgid)
        pnum[hhash] = ii
        ii += 1

ss = 0

ddd = dict()
for key, item in pd.items():
    ddd[item] = ddd.get(item, 0) + 1

print(ddd)
print(len(pd))

# for p in all_particles:
#     hhash = int(p.pdgid) % 1000001
#     print(int(p.pdgid), hhash, pnum[hhash])
#     input()
    # if int(p.pdgid) != ppp[hhash]:
    #     print(int(p.pdgid), ppp[hhash])


# print(pd)
# print([int(p.pdgid) for p in all_particles])

# print(Particle.from_pdgid(PDGID(21212)).name)
#     db = {p.name: p.pdgid for p in all_particles}
#     db.update({p.programmatic_name: p.pdgid for p in all_particles})
#     db["p"] = PDGID(2212)
#     db["n"] = PDGID(2112)
#     db["p~"] = -db["p"]
#     db["n~"] = -db["n"]
#     db.update(
#         H=db["p"],
#         H1=db["p"],
#         He=db["He4"],
#         C=db["C12"],
#         N=db["N14"],
#         O=db["O16"],
#         Ne=db["Ne20"],
#         Ar=db["Ar40"],
#         Xe=db["Xe131"],
#         Pb=db["Pb206"],
#         photon=db["gamma"],
#         proton=db["p"],
#         neutron=db["n"],
#         antiproton=-db["p"],
#         antineutron=-db["n"],
#         pbar=-db["p"],
#         nbar=-db["n"],
#         p_bar=-db["p"],
#         n_bar=-db["n"],
#     )
#     return db
# chromo


nmax = 10000
arr = np.empty(2*nmax + 1, dtype=np.int32)

# arr[5000] = 111
# arr[5001] = 222

# print(arr[-5001], arr[-5000])


icount = 0
for p in all_particles:
    pid = int(p.pdgid)
    if abs(pid) <= nmax:
        arr[pid] = icount
        icount += 1
        
print(f"icount = {icount}")        
# for i in range(-5000, 5000):
#     print(i, arr[i])        


print(arr[111])
print(arr[2212])
print(arr[-2212])


# def pdg2AZ(pdgid):
#     p = PDGID(pdgid)
#     if p.is_nucleus:
#         return p.A, p.Z
#     elif pdgid == 2112:
#         return 1, 0
#     elif pdgid == 2212:
#         return 1, 1
#     return 0, 0

val = 1001172940

print(((val//10) % 1000))
print((val//10000) % 1000)

print(np.int16(32_767))