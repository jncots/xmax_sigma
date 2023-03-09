import MCEq.particlemanager as pman
import MCEq.data as data
import mceq_config as config

config.adv_set["no_mixing"] = True


hdf5_back = data.HDF5Backend()
ics = data.InteractionCrossSections(hdf5_back)
# print(ff)
part = pman.MCEqParticle(
    pdg_id=2212, helicity=0, cs_db=ics, energy_grid=hdf5_back.energy_grid
)
part._calculate_mixing_energy()
# inter = data.Interactions(data.HDF5Backend())
part.set_cs(ics)
print(len(part.inverse_interaction_length()))
