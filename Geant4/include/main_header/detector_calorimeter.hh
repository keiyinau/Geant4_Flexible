#ifndef DETECTOR_CALORIMETER_HH
#define DETECTOR_CALORIMETER_HH

#include "G4ParticleDefinition.hh"
#include "G4VSensitiveDetector.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "run.hh"
#include <map>
#include "SiPMProperties.h"
#include "SiPMAnalogSignal.h"
#include "SiPMSensor.h"
#include <vector>
#include <fstream>  // Required for std::ifstream
#include <iostream>
class Calorimeter : public G4VSensitiveDetector
{
public:
	Calorimeter(G4String);
	~Calorimeter();
	virtual void Initialize(G4HCofThisEvent*); // Add Initialize method
    virtual void EndOfEvent(G4HCofThisEvent*); // For deferred NTuple filling
	void SaveToStepData(G4Step* aStep, G4TouchableHistory* ROhist, G4Track* track);
	void SaveToRoot();
	void ReadOut(G4Step* aStep, G4Track* track);
	void ClearVectorsCounts(); // Clear photon counts and stored data
	struct StepData {
		G4int eventID;
		G4int trackID;
		G4String detectorName;
		G4double scintillatorCount; // Number of scintillators hit by the optical photon
		G4double Hittime;
	};
	std::vector<StepData> CurrentData; // Store exit data for each track
private:
	std::map<G4int, G4double> scintillatorCount; // Map of TrackID to scintillator count
	std::map<G4int, G4double> HitTime; // Map of TrackID to hit time
	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	G4int fHitsCollectionID; // Declare fHitsCollectionID
};

#endif