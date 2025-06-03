#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4ParticleDefinition.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "run.hh"

class MySensentiveDetector : public G4VSensitiveDetector
{
public:
	MySensentiveDetector(G4String);
	~MySensentiveDetector();

	void SaveToDataFile(G4Step* aStep, G4TouchableHistory* ROhist, G4Track* track);
	void ReadOut(G4Step* aStep, G4Track* track);

private:
	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
};

#endif