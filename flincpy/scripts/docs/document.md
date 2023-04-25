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



# Results


## Fluka
<table><tr>
<td> <img src="fig01.png" alt="Muon spectra from Fluka lin" style="width: 400px;" /> </td>
<td> <img src="fig02.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
</tr></table>
<figcaption>

Spectra (histogram) of muons at slant depths $X = (143, 647, 1033)\, g/cm^2$ from air shower initiated by proton with energy $E_0 = 100\text{ GeV}$ and incident zenith angle $\theta = 30^{\circ}$. **Fluka calculations**. "cur" is calculated with currently distributed version - fluka2021.2. "dev" is for the ongoing
development version. Left figure is in linear scale, right figure is in log scale. 

</figcaption>

<table><tr>
<td> <img src="fig03.png" alt="df" style="width: 400px;" /> </td>
</tr></table>
<figcaption>
Ratio of spectra: develop/current 

</figcaption>


## Fluka vs MCEq2D

<table><tr>
<td> <img src="fig04.png" alt="Muon spectra from Fluka lin" style="width: 400px;" /> </td>
<td> <img src="fig05.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
</tr></table>
<figcaption>

Comparison of Fluka(2021) vs MCEq2D energy spectra of muons

</figcaption>

<table><tr>
<td> <img src="fig06.png" alt="df" style="width: 400px;" /> </td>
</tr></table>
<figcaption>
Ratio of spectra: Fluka/MCEq2D

</figcaption>

## Corsika vs MCEq2D

<table><tr>
<td> <img src="fig07.png" alt="Muon spectra from Fluka lin" style="width: 400px;" /> </td>
<td> <img src="fig08.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
</tr></table>
<figcaption>

Comparison of Corsika (Fluka for E < 100 GeV) vs MCEq2D energy spectra of muons

</figcaption>

<table><tr>
<td> <img src="fig09.png" alt="df" style="width: 400px;" /> </td>
</tr></table>
<figcaption>
Ratio of spectra: Corsika/MCEq2D

</figcaption>


# Angular distributions

## Fluka vs MCEq

### X = 143 g/cm2
<table><tr>
<td> <img src="fig10.png" alt="Muon spectra from Fluka lin" style="width: 400px;" /> </td>
<td> <img src="fig11.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
<td> <img src="fig12.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
</tr></table>


### X = 647 g/cm2
<table><tr>
<td> <img src="fig13.png" alt="Muon spectra from Fluka lin" style="width: 400px;" /> </td>
<td> <img src="fig14.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
<td> <img src="fig15.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
</tr></table>

### X = 1033 g/cm2

<table><tr>
<td> <img src="fig16.png" alt="Muon spectra from Fluka lin" style="width: 400px;" /> </td>
<td> <img src="fig17.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
<td> <img src="fig18.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
</tr></table>

---
## Corsika vs MCEq

### X = 143 g/cm2
<table><tr>
<td> <img src="fig19.png" alt="Muon spectra from Fluka lin" style="width: 400px;" /> </td>
<td> <img src="fig20.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
<td> <img src="fig21.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
</tr></table>


### X = 647 g/cm2
<table><tr>
<td> <img src="fig22.png" alt="Muon spectra from Fluka lin" style="width: 400px;" /> </td>
<td> <img src="fig23.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
<td> <img src="fig24.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
</tr></table>

### X = 1033 g/cm2
<table><tr>
<td> <img src="fig25.png" alt="Muon spectra from Fluka lin" style="width: 400px;" /> </td>
<td> <img src="fig26.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
<td> <img src="fig27.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
</tr></table>
---

## Corsika vs MCEq (angular bins as in article)

### X = 143 g/cm2
<table><tr>
<td> <img src="fig28.png" alt="Muon spectra from Fluka lin" style="width: 400px;" /> </td>
<td> <img src="fig29.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
<td> <img src="fig30.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
</tr></table>

### X = 647 g/cm2
<table><tr>
<td> <img src="fig32.png" alt="Muon spectra from Fluka lin" style="width: 400px;" /> </td>
<td> <img src="fig33.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
<td> <img src="fig34.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
</tr></table>

### X = 1033 g/cm2
<table><tr>
<td> <img src="fig35.png" alt="Muon spectra from Fluka lin" style="width: 400px;" /> </td>
<td> <img src="fig36.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
<td> <img src="fig37.png" alt="Muon spectra from Fluka log" style="width: 400px;"/> </td>
</tr></table>
<figcaption>