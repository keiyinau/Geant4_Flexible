// BirksOnsagerSaturation.hh - Updated to match Geant4 11.3.0 signatures
#ifndef BIRKSONSAGERSATURATION_HH
#define BIRKSONSAGERSATURATION_HH

#include "G4EmSaturation.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"

class BirksOnsagerSaturation : public G4EmSaturation {
public:
    BirksOnsagerSaturation(G4int verbose = 0);
    virtual ~BirksOnsagerSaturation();

    // Override the key method to apply BO quenching (add 'const' to match base class)
    virtual G4double VisibleEnergyDeposition(const G4ParticleDefinition* p,
                                             const G4MaterialCutsCouple* couple,
                                             G4double length,
                                             G4double edepTotal,
                                             G4double edepNIEL = 0.0) const override;

    // Optional: Keep this as a convenience method, but remove 'override' since it's not in base
    G4double VisibleEnergyDeposition(const G4Step* step) const;

private:
    // Helper to compute quenching factor L for a given material and dE/dx
    G4double ComputeQuenchingFactor(const G4Material* material, G4double dEdx) const;
};

#endif