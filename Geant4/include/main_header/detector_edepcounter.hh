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
    virtual void Initialize(G4HCofThisEvent*); 
    virtual void EndOfEvent(G4HCofThisEvent*); 
    void SaveToRoot();
    void ClearVectorsCounts(); 

    // Maps grouped ONLY by Detector Name
    std::map<G4String, G4double> edep_per_detector; 
    std::map<G4String, G4double> first_time_per_detector; 
    std::map<G4String, G4int> optical_per_detector; 

private:
    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    G4int fHitsCollectionID; 
};

#endif