#ifndef DETECTOR_TRACKER_HH
#define DETECTOR_TRACKER_HH

#include "G4ParticleDefinition.hh"
#include "G4VSensitiveDetector.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "run.hh"
#include <map>
#include <vector>

class Tracker : public G4VSensitiveDetector
{
public:
	Tracker(G4String);
	~Tracker();
	virtual void Initialize(G4HCofThisEvent*); // Add Initialize method
    virtual void EndOfEvent(G4HCofThisEvent*); // For deferred NTuple filling
	void SaveToStepData(G4Step* aStep, G4TouchableHistory* ROhist, G4Track* track);
	void SaveToRoot();
	void ReadOut(G4Step* aStep, G4Track* track);
	void ClearVectorsCounts(); // Clear photon counts and stored data
	struct StepData {
		G4int eventID;
		G4int trackID;
		G4int stepID;
		G4String particleName;
		G4String ProcessName;
		G4double kineticEnergy;
		G4double AccumatedDistance;
		G4double AccumulatedTime;
		G4double AccumulatedEnergy;

	};
	std::vector<StepData> CurrentData; // Store exit data for each track
private:
    std::map<G4int, G4double> AccumatedDistance_count; // Map of TrackID to photon count
    std::map<G4int, G4double> AccumulatedTime_count; // Map of TrackID to PMT name
	std::map<G4int, G4double> AccumulatedEnergy_count; // Map of TrackID to PMT name
	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	G4int fHitsCollectionID; // Declare fHitsCollectionID
};

#endif