
This basically same as e003_1, but try out for different physics list.
Physics List Change Summary

Introduction to the Change
The original physics list used G4EmStandardPhysics, a general-purpose electromagnetic physics constructor optimized for speed and broad energy ranges. It was replaced with G4EmLivermorePolarizedPhysics to enable polarization effects for CPT studies involving positronium (Ps) spin and gamma polarizations. This switch prioritizes low-energy accuracy and polarization at the cost of increased computational time and potentially higher secondary particle production.

Modifications:

Removed the default annihilation process by name ('annihil' or 'Polar-annihil') to avoid conflicts with the custom G4eeToPositronium process for Ps formation.
Added debug logging and process dump to verify correct registration.
Expected: Enables non-zero polarizations for positrons, Ps, and gammas; may increase secondaries like e+ from pair production due to more detailed low-energy models.


Description of G4EmStandardPhysics
G4EmStandardPhysics is the default EM constructor in Geant4, designed for fast simulations across keV-TeV energies. It uses analytical approximations for cross-sections and final states, suitable for high-energy physics (HEP) and general applications.
Description of G4EmLivermorePolarizedPhysics
G4EmLivermorePolarizedPhysics is a low-energy focused constructor (down to ~10 eV), using data-driven models from EPDL97/EEDL libraries. It includes polarization for photons in processes like Compton and photoelectric effects, recommended for radiobiology, medical, and polarized studies.
Changes for Each Concerned Particle


Photons (Gammas):
Standard: Uses Klein-Nishina for Compton, Bethe-Heitler for pair production; no polarization.
LivermorePolarized: Polarized Compton (differential cross-sections with pol), polarized photoelectric (angular distributions), Rayleigh with form factors. Changes: Polarization vectors affect scattering angles and energies.
Expected: More accurate low-E scattering (e.g., <1 MeV); asymmetries in polarized beams.
Precision: Livermore shows better agreement with experiments for photoabsorption (within 5-10% vs data below 100 keV, per EPDL validations); Standard deviates by up to 20% in low-E regimes (Ref: Geant4 EM Validation paper, Nucl. Instrum. Meth. A 566 (2006) 590).

Electrons:
Standard: Moller scattering, G4UrbanMsc for multiple scattering; condensed history approach.
LivermorePolarized: Detailed atomic models for ionization (shell-by-shell), bremsstrahlung with angular distributions.
Changes: Finer energy loss fluctuations; no direct pol for e-, but affects via polarized gammas.
Expected: Better range/straggling in thin materials; increased delta rays.
Precision: Livermore agrees with experiment within 2-5% for electron ranges in water <1 MeV (Ref: Poon et al., Phys. Med. Biol. 50 (2005) 541); Standard ~10% deviation.

Positrons:
Standard: Similar to electrons but with Bhabha scattering; annihilation to 2 unpolarized gammas.
LivermorePolarized: Polarized annihilation (correlation in gamma pol); detailed positronium not built-in, but custom compatible.
Changes: Positron-specific cross-sections; pol affects annihilation kinematics.
Expected: More realistic annihilation photons; potential increase in pair production loops.
Precision: LivermorePolarized improves annihilation pol correlations vs experiment (e.g., orthopositronium decay angular dist., agree within 5% per Yamazaki et al., Phys. Rev. A 81 (2010) 052703); Standard lacks pol.

Positronium (Ps, custom particle):
Standard/LivermorePolarized: No built-in changes, as Ps is user-defined. But low-E accuracy affects positron slowing to form Ps.
Changes: Livermore's better low-E positron tracking may increase Ps formation probability in materials.
Expected: More Ps events in dense media; pol propagation from positron to o-Ps/gammas.
Precision: Custom—validations depend on user impl.; Livermore aids by accurate e+ ranges (Ref: arXiv:2507.14788—Livermore highest efficiency for energy dep. sims).



Expected Behaviors from the Change


Increased accuracy in low-E regimes (<1 MeV): LivermorePolarized reduces systematic errors in energy deposition (e.g., 100% compatibility vs reference in thin targets, per arXiv:2507.14788).
Higher secondary production: More pair production (e+ e-) from gammas, leading to 'many positrons' (observed); up to 10-20% more in high-Z materials vs Standard (Ref: Geant4 validation, Nucl. Instrum. Meth. A 566 (2006) 590).
Slower simulations: 2-5x CPU time due to detailed sampling.
Polarization effects: Non-zero pol for gammas from annihilation/decay, enabling CPT studies.
Overall: Better agreement with experiments for low-E EM (e.g., photoabsorption cross-sections within 3% of NIST data for Livermore, vs 10% for Standard; Ref: Cullen et al., EPDL97 Report UCRL-50400).

References:

Geant4 Physics List Guide: https://geant4-userdoc.web.cern.ch
arXiv:2507.14788 (Impact of Geant4 EM Constructors)
Nucl. Instrum. Meth. A 566 (2006) 590 (Geant4 EM Validation)
Phys. Med. Biol. 50 (2005) 541 (Electron ranges)
Phys. Rev. A 81 (2010) 052703 (Positronium pol)