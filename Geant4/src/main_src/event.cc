#include "event.hh"

MyEventAction::MyEventAction(MyRunAction*)
{

}

MyEventAction::~MyEventAction()
{}

void MyEventAction::BeginOfEventAction(const G4Event* aEvent)
{

	G4cout << ">> Begin of Event:" << aEvent->GetEventID() << G4endl;
    positronCreators.clear();
    psPositions.clear(); psMomenta.clear(); psPols.clear(); psTypes.clear(); psParents.clear(); psDestroyTimes.clear();
    gammaPositions.clear(); gammaMomenta.clear(); gammaPols.clear(); gammaEnergies.clear(); gammaTypes.clear(); gammaParents.clear(); gammaFirstDets.clear();
    positronPositions.clear(); positronMomenta.clear(); positronPols.clear();
    primaryDecayTime = 0.0;

}
void MyEventAction::SetPrimaryDecayTime(G4double t) {
    if (primaryDecayTime == 0.0) {  // Set only once (primary)
        primaryDecayTime = t;
    }
}
void MyEventAction::EndOfEventAction(const G4Event* aEvent)
{
G4int evt = aEvent->GetEventID();
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    for (const auto& entry : positronPositions) {
        G4int trk = entry.first;
        man->FillNtupleIColumn(6, 0, evt);
        man->FillNtupleIColumn(6, 1, trk);
        man->FillNtupleDColumn(6, 2, positronPositions[trk].x() / mm);
        man->FillNtupleDColumn(6, 3, positronPositions[trk].y() / mm);
        man->FillNtupleDColumn(6, 4, positronPositions[trk].z() / mm);
        man->FillNtupleDColumn(6, 5, positronMomenta[trk].x() / MeV);
        man->FillNtupleDColumn(6, 6, positronMomenta[trk].y() / MeV);
        man->FillNtupleDColumn(6, 7, positronMomenta[trk].z() / MeV);
        man->FillNtupleDColumn(6, 8, positronPols[trk].x());
        man->FillNtupleDColumn(6, 9, positronPols[trk].y());
        man->FillNtupleDColumn(6, 10, positronPols[trk].z());
        man->FillNtupleSColumn(6, 11, positronCreators[trk]);
        man->FillNtupleDColumn(6, 12, positronTimes[trk] / ns);
        man->AddNtupleRow(6);
        positronTime = positronTimes[trk]- primaryDecayTime ;
    }
    for (const auto& entry : psPositions) {
        G4int trk = entry.first;
        man->FillNtupleIColumn(4, 0, evt);
        man->FillNtupleIColumn(4, 1, trk);
        man->FillNtupleIColumn(4, 2, psParents[trk]);
        man->FillNtupleSColumn(4, 3, psTypes[trk]);
        man->FillNtupleDColumn(4, 4, psPositions[trk].x() / mm);
        man->FillNtupleDColumn(4, 5, psPositions[trk].y() / mm);
        man->FillNtupleDColumn(4, 6, psPositions[trk].z() / mm);
        man->FillNtupleDColumn(4, 7, psMomenta[trk].x() / MeV);
        man->FillNtupleDColumn(4, 8, psMomenta[trk].y() / MeV);
        man->FillNtupleDColumn(4, 9, psMomenta[trk].z() / MeV);
        man->FillNtupleDColumn(4, 10, psPols[trk].x());
        man->FillNtupleDColumn(4, 11, psPols[trk].y());
        man->FillNtupleDColumn(4, 12, psPols[trk].z());
        G4double lifetime = (psDestroyTimes[trk]-primaryDecayTime - positronTime) / ns;
        man->FillNtupleDColumn(4, 13, lifetime);
        man->AddNtupleRow(4);
    }

    for (const auto& entry : gammaPositions) {
        G4int trk = entry.first;
        G4String firstDet = (gammaFirstDets.count(trk) > 0) ? gammaFirstDets[trk] : "None";
        man->FillNtupleIColumn(5, 0, evt);
        man->FillNtupleIColumn(5, 1, trk);
        man->FillNtupleIColumn(5, 2, gammaParents[trk]);
        man->FillNtupleSColumn(5, 3, gammaTypes[trk]);
        man->FillNtupleDColumn(5, 4, gammaEnergies[trk] / MeV);
        man->FillNtupleDColumn(5, 5, gammaPositions[trk].x() / mm);
        man->FillNtupleDColumn(5, 6, gammaPositions[trk].y() / mm);
        man->FillNtupleDColumn(5, 7, gammaPositions[trk].z() / mm);
        man->FillNtupleDColumn(5, 8, gammaMomenta[trk].x() / MeV);
        man->FillNtupleDColumn(5, 9, gammaMomenta[trk].y() / MeV);
        man->FillNtupleDColumn(5, 10, gammaMomenta[trk].z() / MeV);
        man->FillNtupleDColumn(5, 11, gammaPols[trk].x());
        man->FillNtupleDColumn(5, 12, gammaPols[trk].y());
        man->FillNtupleDColumn(5, 13, gammaPols[trk].z());
        man->FillNtupleSColumn(5, 14, firstDet);
        man->AddNtupleRow(5);
    }


    for (const auto& entry : gammaEdeps) {
        G4int trk = entry.first;
        for (const auto& edep : entry.second) {
            man->FillNtupleIColumn(7, 0, evt);
            man->FillNtupleIColumn(7, 1, trk);
            man->FillNtupleDColumn(7, 2, edep.deltaE / MeV);
            man->FillNtupleSColumn(7, 3, edep.detName);
            man->FillNtupleDColumn(7, 4, (edep.time - primaryDecayTime) / ns);
            man->AddNtupleRow(7);
        }
    }
}

void MyEventAction::AddPsTruth(G4int trackID, G4int parentID, G4String type, G4ThreeVector pos, G4ThreeVector mom, G4ThreeVector pol) {
    psPositions[trackID] = pos;
    psMomenta[trackID] = mom;
    psPols[trackID] = pol;
    psTypes[trackID] = type;
    psParents[trackID] = parentID;
}

void MyEventAction::AddGammaTruth(G4int trackID, G4int parentID, G4String type, G4double energy, G4ThreeVector pos, G4ThreeVector mom, G4ThreeVector pol) {
    gammaPositions[trackID] = pos;
    gammaMomenta[trackID] = mom;
    gammaPols[trackID] = pol;
    gammaEnergies[trackID] = energy;
    gammaTypes[trackID] = type;
    gammaParents[trackID] = parentID;
}

void MyEventAction::AddPositronTruth(G4int trackID, G4ThreeVector pos, G4ThreeVector mom, G4ThreeVector pol,G4String creatorProcess, G4double time) {
    positronPositions[trackID] = pos;
    positronMomenta[trackID] = mom;
    positronPols[trackID] = pol;
    positronCreators[trackID] = creatorProcess;
    positronTimes[trackID] = time;
}

void MyEventAction::SetGammaFirstDetector(G4int trackID, G4String detName) {
    if (gammaFirstDets.count(trackID) == 0) { // Only set first
        gammaFirstDets[trackID] = detName;
    }
}
void MyEventAction::AddGammaEdep(G4int trackID, G4double deltaE, G4String detName, G4double time) {
    if (gammaPositions.count(trackID) > 0) {  // Only for TruthGammas
        GammaEdep edepInfo = {deltaE, detName, time};
        gammaEdeps[trackID].push_back(edepInfo);
    }
}
void MyEventAction::AddPsDestroyTime(G4int trackID, G4double time) {
    psDestroyTimes[trackID] = time;
}