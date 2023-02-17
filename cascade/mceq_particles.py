from MCEq.particlemanager import MCEqParticle, ParticleManager
from MCEq.data import InteractionCrossSections, HDF5Backend


hdf5_backend = HDF5Backend()
interaction_cs = InteractionCrossSections(hdf5_backend)


print(interaction_cs.get_cs(2212))



# pman = ParticleManager()











