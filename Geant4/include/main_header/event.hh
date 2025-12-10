#ifndef EVENT_HH
#define EVENT_HH

#include <map>
#include <string>

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"

#include "run.hh"

class MyEventAction : public G4UserEventAction
{
public:
	MyEventAction(MyRunAction*);
    ~MyEventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    // New methods for truth collection
    void AddPsTruth(G4int trackID, G4int parentID, G4String type, G4ThreeVector pos, G4ThreeVector mom, G4ThreeVector pol);
    void AddGammaTruth(G4int trackID, G4int parentID, G4String type, G4double energy, G4ThreeVector pos, G4ThreeVector mom, G4ThreeVector pol);
    void AddPositronTruth(G4int trackID, G4ThreeVector pos, G4ThreeVector mom, G4ThreeVector pol);
    void SetGammaFirstDetector(G4int trackID, G4String detName);
    G4bool HasPs(G4int id) { return psPositions.count(id) > 0; }
private:
	G4double fEdep;
	std::map<G4int, G4ThreeVector> psPositions, psMomenta, psPols;
    std::map<G4int, G4String> psTypes;
    std::map<G4int, G4int> psParents;

    std::map<G4int, G4ThreeVector> gammaPositions, gammaMomenta, gammaPols;
    std::map<G4int, G4double> gammaEnergies;
    std::map<G4int, G4String> gammaTypes;
    std::map<G4int, G4int> gammaParents;
    std::map<G4int, G4String> gammaFirstDets; // First detector per gamma trackID

    std::map<G4int, G4ThreeVector> positronPositions, positronMomenta, positronPols;
};
#endif