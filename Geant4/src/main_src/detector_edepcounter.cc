#include "detector_edepcounter.hh"
#include "event.hh"

Detect_edep::Detect_edep(G4String name) : G4VSensitiveDetector(name), fHitsCollectionID(-1)
{
    ClearVectorsCounts(); 
    collectionName.insert("EdepCollection");
}

Detect_edep::~Detect_edep() {}

void Detect_edep::Initialize(G4HCofThisEvent* hce)
{
    ClearVectorsCounts(); // Clear everything at the start of the event
    
    if (fHitsCollectionID < 0) {
        fHitsCollectionID = GetCollectionID(0);
    }
    G4VHitsCollection* hc = new G4VHitsCollection(SensitiveDetectorName, collectionName[0]);
    hce->AddHitsCollection(fHitsCollectionID, hc);
}

void Detect_edep::EndOfEvent(G4HCofThisEvent*)
{
    SaveToRoot();
    ClearVectorsCounts(); 
}

G4bool Detect_edep::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    G4Track* track = aStep->GetTrack();
    G4String detector_Name = track->GetTouchable()->GetVolume()->GetName();
    G4String particle = track->GetParticleDefinition()->GetParticleName();
    G4double edep_step = aStep->GetTotalEnergyDeposit();

    // Handle Optical Photons separately
    if (particle == "opticalphoton") {
        if (track->GetCurrentStepNumber() == 1) {
            optical_per_detector[detector_Name]++;
        }
        return true;
    }

    // Handle standard particles (electrons, gammas, etc.) depositing energy
    if (edep_step > 0.) {
        
        // 1. Accumulate total energy into this detector
        edep_per_detector[detector_Name] += edep_step;

        // 2. Record the earliest time of interaction in this detector
        G4double time = aStep->GetPreStepPoint()->GetGlobalTime();
        
        // If this detector hasn't been hit yet, OR if this hit is earlier than the saved one
        if (first_time_per_detector.find(detector_Name) == first_time_per_detector.end() || time < first_time_per_detector[detector_Name]) {
            first_time_per_detector[detector_Name] = time;
        }
    }

    return true;
}

void Detect_edep::SaveToRoot()
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    const MyEventAction* eventAction = static_cast<const MyEventAction*>(G4RunManager::GetRunManager()->GetUserEventAction());
    
    G4double primaryTime = 0.0;
    if (!eventAction) {
        G4cerr << "Warning: No EventAction found — using absolute times" << G4endl;
    } else {
        primaryTime = eventAction->GetPrimaryDecayTime();
    }

    // Loop through each detector that had energy deposited
    for (const auto& pair : edep_per_detector) {
        G4String detName = pair.first;
        G4double edep_acc = pair.second;
        
        // Only save if energy was actually deposited
        if (edep_acc > 0. * eV) {
            analysisManager->FillNtupleIColumn(1, 0, evt); // eventID
            analysisManager->FillNtupleSColumn(1, 1, detName); // detectorName
            
            // The TOTAL accumulated edep for the whole detector
            analysisManager->FillNtupleDColumn(1, 2, edep_acc / MeV); 
            
            // Calculate the relative time using the earliest hit in the detector
            G4double firstTime = first_time_per_detector[detName];
            G4double rel_time = (firstTime - primaryTime) / ns;
            
            analysisManager->FillNtupleDColumn(1, 3, rel_time); 
            analysisManager->FillNtupleIColumn(1, 4, optical_per_detector[detName]);
            analysisManager->AddNtupleRow(1);
        }
    }
}

void Detect_edep::ClearVectorsCounts()
{
    edep_per_detector.clear();
    first_time_per_detector.clear();
    optical_per_detector.clear();
}