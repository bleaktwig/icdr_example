# Data release of the IceCube/DeepCore 9 year golden event oscillation analysis
This sample contains the data used in the [Measurement of atmospheric neutrino mixing with improved IceCube DeepCore calibration and data processing](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.108.012014), published in Phys. Rev. D 108, 012014 (2023). Please refer to the publication for a detailed explanation of the sample and analysis.

## Files
Along with this readme, you'll find six `.csv` files and one ipython notebook:
* `data.csv`: count of data events per analysis bin.
* `mc_nu*.csv`: event-by-event information of MC neutrinos, split between neutral current (NC) and particle type + charged current (CC).
* `mc_mu.csv`: per-bin information from the MC muon background.
* `example.ipynb`: ipython notebook presenting a simple example fit to the data.

## Binning
Reconstructed variables are provided in the analysis binning, such that values are chosen so that events fall into the correct bins.
* Reconstructed energy or `reco_energy` in GeV: [6.31, 8.45862141, 11.33887101, 15.19987592, 20.37559363, 27.3136977, 36.61429921, 49.08185342, 65.79474104, 88.19854278, 158.49].
* Reconstructed cosine of the particle's zenith angle or `reco_coszen`: [-1., -0.89, -0.78, -0.67, -0.56, -0.45, -0.34, -0.23, -0.12, -0.01,  0.1]. The convention is such that -1 corresponds to straight "upgoing" trajectories (earth crossing), +1 to straight "downgoing" (directly from the sky above), and 0 is horizontal.
* Reconstructed particle ID or `pid`: [0.55, 0.75, 1.]. Represents the likelihood that an event is a cascade or a track: 0 corresponds to a cascade, and 1 to a track to the best of our knowledge.

---
## data.csv
Counts of observed data in the analysis binning. Each row corresponds to one bin, such that
* `count`: count of data events in the corresponding bin.
* `pid`: reconstructed particle ID in analysis binning.
* `reco_coszen`: reconstructed cosine(zenith) in analysis binning.
* `reco_energy`: reconstructed energy (GeV) in analysis binning.

---
## mc_nu*.csv
MC neutrinos with one event per row. The four files included are
* `mc_nu_nc.csv`: neutral-current (NC) neutrino events, all flavors;
* `mc_nue_cc.csv`: charged-current (CC) electron neutrino events;
* `mc_numu_cc.csv`: CC muon neutrino events; and
* `mc_nutau_cc.csv`: CC tau neutrino events.

Each row correspond to one event. The reconstructed variables follow the same convention as `data.csv`, where
* `reco_energy` is the reconstructed energy (GeV) in analysis binning;
* `reco_coszen` is the reconstructed cosine(zenith) in analysis binning; and
* `pid` is the reconstructed particle ID in analysis binning.

Then, the true information is
* `pdg` stands for the PDG code of the neutrino interacting in the detector (see the [PDG summary tables](https://pdg.lbl.gov/2024/tables/contents_tables.html));
* `true_energy` is the simulated energy (GeV) of the neutrino; and
* `true_coszen` is the simulated cosine(zenith) of the neutrino.

Due to the way neutrinos are simulated in IceCube, each MC event has an associated weight, encoded by
* `weight`, in (GeV cm^2 sr). This should be multiplied by the atmospheric flux in (cm^-2 sr^-1 s^-1) and by livetime (s). For the sample's livetime, please see the analysis publication.

There are a few interaction types for the simulated neutrinos, which are encoded into
* `type`, where,
    * 0 stands for any NC event;
    * 1 are CC quasi-elastic (CCQE) events;
    * 2 are CC resonant (CCRES) events;
    * 3 are CC deep inelastic scattering (DIS) events;
    * 4 are CC electromagnetic interaction (CCEM); and
    * -1 stands for other types of interactions. These are very rare (< 0.04% of the total weighted events).

For non-DIS CC events, we encode the interaction uncertainty parameters through two linear fits,
* `MaCCQE` for the axial mass CC quasi-elastic events; and
* `MaCCRES` for the axial mass CC resonant events.

Then, for CC DIS events, we encode the experimental-like kinematic variables. All of these neglect Fermi momentum / off-shellness.
* `Q2` (GeV^2) is the momentum transfer Q^2;
* `W` (GeV) is the hadronic invariant mass;
* `x` is the Bjorken scaling variable; and
* `y` the inelasticity.

Afterwards, there are 10 columns which correspond to the ultrasurfaces gradients -- 2 per detector systematic. These parameters are used to update the weight in each event: if `d_syst` is the difference between the detector systematic's best fit value and the nominal value, `grad_syst` and `grad2_syst` the gradient parameters, and `W` the weight to be updated, then

        W = W * (1 + grad_syst * d_syst + grad2_syst * d_syst^2),

for each detector systematic. The nominal detector systematic parameters used for this study are
* DOM efficiency: 100%.
* Hole ice p0: 0.101569.
* Hole ice p1: -0.049344.
* Bulk ice absorption: 1.
* Bulk ice scattering: 1.05.

For a practical example on how to reweight events, please see the `ultrasurfaces_reweight()` function in the attached notebook. For a detailed explanation of the method, please see [Treating detector systematics via a likelihood free inference method](https://iopscience.iop.org/article/10.1088/1748-0221/18/10/P10019).

Among these columns,
* `grad_dom_eff` and `grad2_dom_eff` are the gradient parameters for the DOM efficiency;
* `grad_hi_p0` and `grad2_hi_p0` are those for the hole ice  p0;
* `grad_hi_p1` and `grad2_hi_p1` for the hole ice p1;
* `grad_bi_abs` and `grad2_bi_abs` for the bulk ice absorption; and
* `grad_bi_scat` and `grad2_bi_scat` for the bulk ice scattering.

---
## mc_mu.csv
Counts and absolute uncertainties for the atmospheric muon background in the analysis binning. Due to the high efficiency of the event selection, the produced histogram is sparsely populated, which is overcome by applying a variable bandwidth kernel density estimator (KDE).
* `count`: count of muons in the corresponding bin.
* `abs_uncertainty`: absolute uncertainty on count, accounting for statistical and shape uncertainty.
* `pid`: reconstructed particle ID in analysis binning.
* `reco_coszen`: reconstructed cosine(zenith) in analysis binning.
* `reco_energy`: reconstructed energy (GeV) in analysis binning.
