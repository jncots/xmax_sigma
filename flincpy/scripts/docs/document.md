# Parameters of a shower
- E = 100 GeV
- Theta = 30 deg

```
Hi Anatoli

Attached please find a text file with 9 angular distributions.

They are for the 3 muon energy intervals, 1-1.3, 2-2.5, 4-5 GeV,
at 143 g.cm^2,  647 g.cm^2, and 1033 g/cm^2 dpeth, for 100 GeV protons
inpinging at 30 degress on the atmosphere.

The numbers for each distributions are:

Omega_i-1   Omega_i    d2N/dOmega/dE     err
   sr          sr      1/sr/GeV/prim      %

You can obviously get the angles out of the Omega solid angles,
from:

  Omega_i = 2 pi [1-cos(theta_i)]

The 9 distributions have name names like Ej-<depth>, where
j=1 is for the 1-1.3 GeV interval, J=2 for the 2-2.5 one, and j=3 for
the 4-5 GeV one.

The statistics is for 43500000 primaries.

Let me know if the numbers are fishy, I did it somewhat in a hurry, I
cannot exclude I made some mistakes.
```

# Filtered data
- Pdg:
    - [-13, 13] muons
    - [-12, 12] nu_e
    - [-14, 14] nu_mu

- Xdepth:
    - 143 g/cm2
    - 647 g/cm2
    - 1033 g/cm2

- Angular dist at :
    - 1-1.3 GeV
    - 2-2.5 GeV
    - 4-5 GeV 


# Models

## Fluka
- Fluka 2021
- Fluka develop

## Corsika (new)
- Fluka 2021

## MCEq (2D)
- Below E_th = 80 GeV
    - UrQMD
    - DPMJet-III 19.1

- Above:
Sibyll-2.3d, EPOS-LHC, DPMJet-III 19.1  
