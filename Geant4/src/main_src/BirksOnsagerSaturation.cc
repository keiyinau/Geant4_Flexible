// BirksOnsagerSaturation.cc - Updated to match Geant4 11.3.0
#include "BirksOnsagerSaturation.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonisParamMat.hh"
#include <cmath>  // for exp

BirksOnsagerSaturation::BirksOnsagerSaturation(G4int verbose) : G4EmSaturation(verbose) {}

BirksOnsagerSaturation::~BirksOnsagerSaturation() {}

G4double BirksOnsagerSaturation::VisibleEnergyDeposition(const G4ParticleDefinition* p,
                                                         const G4MaterialCutsCouple* couple,
                                                         G4double length,
                                                         G4double edepTotal,
                                                         G4double edepNIEL) const {
    if (length <= 0.0 || edepTotal <= 0.0) return edepTotal - edepNIEL;

    // Compute dE/dx (MeV/mm, assuming length in mm)
    G4double dEdx = (edepTotal - edepNIEL) / (length / mm);  // Normalize to MeV/mm

    const G4Material* material = couple->GetMaterial();

    // Compute full BO quenching factor
    G4double L = ComputeQuenchingFactor(material, dEdx);

    // Return quenched visible energy
    return L * (edepTotal - edepNIEL);
}

// Optional convenience method (not overriding base)
G4double BirksOnsagerSaturation::VisibleEnergyDeposition(const G4Step* step) const {
    if (!step) return 0.0;

    G4double edepTotal = step->GetTotalEnergyDeposit();
    G4double length = step->GetStepLength();
    if (length <= 0.0 || edepTotal <= 0.0) return edepTotal;

    // NIEL (non-ionizing) is not directly available in step; assume 0 for simplicity or compute if needed
    G4double edepNIEL = 0.0;  // Customize if you have NIEL modeling

    // dE/dx in MeV/mm
    G4double dEdx = (edepTotal - edepNIEL) / (length / mm);

    const G4Material* material = step->GetPreStepPoint()->GetMaterial();

    // Compute full BO quenching factor
    G4double L = ComputeQuenchingFactor(material, dEdx);

    return L * (edepTotal - edepNIEL);
}

G4double BirksOnsagerSaturation::ComputeQuenchingFactor(const G4Material* material, G4double dEdx) const {
    // Get Birks constant (mm/MeV)
    G4double kB = material->GetIonisation()->GetBirksConstant() / (mm / MeV);

    // Defaults (paper values)
    G4double etaH    = 0.0;
    G4double eta_eh  = 0.0;
    G4double k_O     = 0.29 * mm / MeV;   // ← paper value, mm/MeV

    // Read from material properties table
    G4MaterialPropertiesTable* mpt = material->GetMaterialPropertiesTable();
    if (mpt) {
        if (mpt->ConstPropertyExists("BIRKS_ETA_H")) {
            etaH = mpt->GetConstProperty("BIRKS_ETA_H");
        }
        if (mpt->ConstPropertyExists("ONSAGER_ETA_EH")) {
            eta_eh = mpt->GetConstProperty("ONSAGER_ETA_EH");
        }
        if (mpt->ConstPropertyExists("ONSAGER_K_O")) {
            k_O = mpt->GetConstProperty("ONSAGER_K_O") * (mm / MeV);
        }
    }

    // Birks term (Generalised Birks) — unchanged
    G4double birksTerm = (1.0 - etaH) / (1.0 + kB * (1.0 - etaH) * dEdx) + etaH;

    // Onsager term — exactly as in the paper
    G4double onsagerTerm = 1.0 - eta_eh * std::exp(-k_O * dEdx);

    return onsagerTerm * birksTerm;
}