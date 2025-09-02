#ifndef DETECTOR_IDEAL_HH
#define DETECTOR_IDEAL_HH

#include "G4ParticleDefinition.hh"
#include "G4VSensitiveDetector.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "run.hh"
#include <map>
#include <vector>

class Detect_reference : public G4VSensitiveDetector
{
public:
	Detect_reference(G4String);
	~Detect_reference();
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
		G4int parentID; // Parent ID of the track
		G4String detectorName;
		G4String particleName;
		G4String creatorProcessName; // Process that created the track
		G4String ProcessName;
		G4double kineticEnergy;
		G4double x_distance;
		G4double y_distance;
		G4double z_distance;
		G4double x_momentum;
		G4double y_momentum;
		G4double z_momentum;
	};
	std::vector<StepData> CurrentData; // Store exit data for each track
private:
	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	G4int fHitsCollectionID; // Declare fHitsCollectionID
};

#endif