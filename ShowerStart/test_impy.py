import numpy as np
import impy
from impy.constants import TeV
from impy.kinematics import EventFrame

# Define the parameters of the collisions
event_kinematics = impy.kinematics.FixedTarget(13 * TeV, "proton", "proton")
# Create an instance of an event generator
generator = impy.models.Sibyll23d(event_kinematics)
generator._frame = EventFrame.FIXED_TARGET

nevents = 0
average_pt = 0

# Generate 10000 events
for event in generator(3):
    # Filter event
    event = event.final_state()
    print(f"pid = {event.pid}, energy = {event.en}\n")
    # do something with event.pid, event.eta, event.en, event.pt, etc.
    # these variables are numpy arrays, that can be histogrammed or counted like
    pt = event.pt[np.abs(event.pid) == 211]
    # The list could be empty
    if len(pt) > 0:
        nevents += 1
        average_pt += np.mean(pt)

average_pt = average_pt / nevents
print("Average pT for charged pions {0:4.3f}".format(average_pt))



