#ifndef DETECTOR_EDEPCOUNTER_HH
#define DETECTOR_EDEPCOUNTER_HH

#include "G4ParticleDefinition.hh"
#include "G4VSensitiveDetector.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "run.hh"
#include <map>
#include <vector>

class Detect_edep : public G4VSensitiveDetector
{
public:
    Detect_edep(G4String);
    ~Detect_edep();
    virtual void Initialize(G4HCofThisEvent*); // Add Initialize method
    virtual void EndOfEvent(G4HCofThisEvent*); // For deferred NTuple filling
    void SaveToStepData(G4Step* aStep, G4TouchableHistory* ROhist, G4Track* track);
    void SaveToRoot();
    void ReadOut(G4Step* aStep, G4Track* track);
    void ClearVectorsCounts(); // Clear photon counts and stored data
    struct StepData {
        G4int eventID;
        G4String detectorName;
        G4double edep_accumulated;
    };
    std::map<G4String, G4double> edep_per_detector; // Map to accumulate edep per detector name
    std::map<G4String, G4double> first_time_per_detector; // New: Map for earliest time per detector (in seconds)
    std::vector<StepData> CurrentData; // Store exit data for each track (if needed)
private:
    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    G4int fHitsCollectionID; // Declare fHitsCollectionID
};

#endif