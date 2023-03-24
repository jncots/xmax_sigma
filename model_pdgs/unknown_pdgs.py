from MCEq.core import MCEqRun
import crflux.models as pm
import platform
from tqdm import tqdm

def unique_pdg_list(pdgs):
    res = list(set(pdgs))
    res.sort(key=lambda x: (abs(x), x > 0))
    return res


class MCEqPdgsList:
    def __init__(self):
        self.mceq_run = MCEqRun(
        #provide the string of the interaction model
        interaction_model="DPMJET-III-19.1",
        #primary cosmic ray flux model
        primary_model = (pm.HillasGaisser2012, "H3a"),
        # Zenith angle in degrees. 0=vertical, 90=horizontal
        theta_deg=0.0
        )
    def __call__(self):
        return unique_pdg_list([p.pdg_id[0] for p in self.mceq_run.pman.all_particles])

def unknown_pdgs_list(event_generator, known_pdgs, number_of_events = 10000):
    all_pdgs = set()
    for event in tqdm(event_generator(number_of_events), total = number_of_events):
        all_pdgs.update(event.final_state().pid) 
    return unique_pdg_list(all_pdgs.difference(known_pdgs))


def get_all_models(skip=None):
    from chromo import models
    from chromo.common import MCRun

    if skip is None:
        skip = []
    if platform.system() == "Windows":
        skip = list(skip) + [models.UrQMD34, models.Pythia8]

    result = []
    for key in dir(models):
        obj = getattr(models, key)
        if skip and obj in skip:
            continue
        try:
            if issubclass(obj, MCRun):  # fails if obj is not a class
                result.append(obj)
        except TypeError:
            pass

    return result

    