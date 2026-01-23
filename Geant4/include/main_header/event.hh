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

    void AddPsTruth(G4int trackID, G4int parentID, G4String type, G4ThreeVector pos, G4ThreeVector mom, G4ThreeVector pol,G4double createTime);
    void AddGammaTruth(G4int trackID, G4int parentID, G4String type, G4double energy, G4ThreeVector pos, G4ThreeVector mom, G4ThreeVector pol);
    void AddPositronTruth(G4int trackID, G4ThreeVector pos, G4ThreeVector mom, G4ThreeVector pol, G4String creatorProcess, G4double time);
    void SetGammaFirstDetector(G4int trackID, G4String detName);
    G4bool HasPs(G4int id) { return psPositions.count(id) > 0; }
    struct GammaEdep {
    G4double deltaE;      // MeV
    G4String detName;     // ns (global time)};
    G4double time;        // ns (global time)
    };
    void AddGammaEdep(G4int trackID, G4double deltaE, G4String detName, G4double time);
    std::map<G4int, std::vector<GammaEdep>> gammaEdeps; // trackID -> list of energy depositions
    void SetPrimaryDecayTime(G4double t);
    G4double GetPrimaryDecayTime() const { return primaryDecayTime; }

    void AddPsDestroyTime(G4int trackID, G4double time);
private:
	G4double fEdep;
    G4double primaryDecayTime;
    G4double positronTime;
	std::map<G4int, G4ThreeVector> psPositions, psMomenta, psPols;
    std::map<G4int, G4String> psTypes;
    std::map<G4int, G4int> psParents;
    std::map<G4int, G4double> pscreationtime, pslifetimes;
    std::map<G4int, G4double> psDestroyTimes;
    std::map<G4int, G4double> psCreateTimes;  
    std::map<G4int, G4ThreeVector> gammaPositions, gammaMomenta, gammaPols;
    std::map<G4int, G4double> gammaEnergies;
    std::map<G4int, G4String> gammaTypes;
    std::map<G4int, G4int> gammaParents;
    std::map<G4int, G4String> gammaFirstDets; // First detector per gamma trackID
    
    std::map<G4int, G4ThreeVector> positronPositions, positronMomenta, positronPols;
    std::map<G4int, G4String> positronCreators;
    std::map<G4int, G4double> positronTimes;

};
#endif