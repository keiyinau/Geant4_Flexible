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
#include <fstream> // Required for std::ifstream
#include <iostream>
#include <TGraph.h>
#include <TCanvas.h>
#include <float.h>  // For DBL_MAX (alternatively, #include <limits> and use std::numeric_limits<double>::max())

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
    struct LoadData {
        G4int eventID;
        G4String detectorName;  // Renamed from SiPMName for clarity
        G4double Area;
        G4int RealPhotonCount;
        G4int PEsCount;
        G4int NoisePEsCount;
        G4double Time_Of_Triggering;
    };
    std::vector<LoadData> CurrentData; // Store exit data for each track
    G4String detectorname; // Store exit data for each track (may be removable if unused)
    sipm::SiPMProperties myProperties;
    sipm::SiPMSensor mySensor;
    sipm::SiPMAnalogSignal mySignal;
    void PlotWaveform(const sipm::SiPMAnalogSignal& signal, const G4String& det_name = "");  // Added param
    double signalLength, SampleTime, DarkCountRate, RiseTime, FallTimeFast, RecoveryTime, Dcr, Xt, Ap, pitch, size, gain;
    int nCells;
    double gatewidth, threshold;
    bool isGraph, isDCR, isXT, isAP;

private:
    std::map<G4String, std::vector<double>> photonTimes_per_detector;  // New: per detector
    std::map<G4String, std::vector<double>> photonWavelengths_per_detector;  // New: per detector
    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    G4int fHitsCollectionID; // Declare fHitsCollectionID
};
#endif