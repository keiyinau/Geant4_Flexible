#include "detector_edepcounter.hh"

Detect_edep::Detect_edep(G4String name) : G4VSensitiveDetector(name), fHitsCollectionID(-1)
{
    ClearVectorsCounts(); // Initialize the vectors to store accumulated data
    collectionName.insert("EdepCollection");
}

Detect_edep::~Detect_edep()
{}

void Detect_edep::Initialize(G4HCofThisEvent* hce)
{
    edep_per_detector.clear(); // Clear the map at the start of each event
    first_time_per_detector.clear(); // New: Clear the time map
    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    //G4cout << "MySensitiveDetector::Initialize called for Event=" << eventID << G4endl;
    if (fHitsCollectionID < 0) {
        fHitsCollectionID = GetCollectionID(0);
    }
    G4VHitsCollection* hc = new G4VHitsCollection(SensitiveDetectorName, collectionName[0]);
    hce->AddHitsCollection(fHitsCollectionID, hc);
}

void Detect_edep::EndOfEvent(G4HCofThisEvent*)
{
    SaveToRoot();
    ClearVectorsCounts(); // Clear the accumulated counts at the end of each event
}

G4bool Detect_edep::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
    G4Track* track = aStep->GetTrack();
    G4String detector_Name = track->GetTouchable()->GetVolume()->GetName();
    G4String particle = track->GetParticleDefinition()->GetParticleName();
    G4double edep_step = aStep->GetTotalEnergyDeposit();

    if (particle != "opticalphoton" && edep_step > 0.) { // Skip optical photons and zero-edep steps
        edep_per_detector[detector_Name] += edep_step;

        // New: Record the earliest global time for the first interaction (min time of depositing steps)
        G4double time = aStep->GetPreStepPoint()->GetGlobalTime();
        auto it = first_time_per_detector.find(detector_Name);
        if (it == first_time_per_detector.end() || time < it->second) {
            first_time_per_detector[detector_Name] = time;
        }
    }
    // Optionally call ReadOut(aStep, track) for debugging
    return true;
}

// This function stores information to an Ntuple then it can be saved in run.cc
void Detect_edep::SaveToStepData(G4Step* aStep, G4TouchableHistory* ROhist, G4Track* track)
{
    // This function is no longer needed for per-step saving; accumulation happens in ProcessHits
}

void Detect_edep::SaveToRoot()
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();


    G4double min_time = DBL_MAX; // Use a large initial value
    for (const auto& pair : first_time_per_detector) {
        if (pair.second < min_time) {
            min_time = pair.second;
        }
    }
    // If no times recorded (no deposits), set min_time to 0
    if (min_time == DBL_MAX) {
        min_time = 0.;
    }

    for (const auto& pair : edep_per_detector) {
        if (pair.second >= 500. * eV) {
            analysisManager->FillNtupleIColumn(1, 0, evt); // eventID
            analysisManager->FillNtupleSColumn(1, 1, pair.first); // detectorName
            analysisManager->FillNtupleDColumn(1, 2, pair.second / MeV); // edep_accumulated
            // New: Fill the first time (in ns; adjust unit if needed)
            G4double rel_time = (first_time_per_detector[pair.first] - min_time) / ns;
            analysisManager->FillNtupleDColumn(1, 3, rel_time);
            analysisManager->AddNtupleRow(1);
        }
    }
}

// Output Information just touch the detector
void Detect_edep::ReadOut(G4Step* step, G4Track* track) {

    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    G4int trackID = track->GetTrackID();
    G4int stepID = track->GetCurrentStepNumber();
    G4String particle_name = track->GetDefinition()->GetParticleName();
    G4String creator_process_name = "NULL";
    G4String physVol_name = track->GetTouchable()->GetVolume()->GetName();
    G4ThreeVector postDetectorPosition = track->GetTouchable()->GetVolume()->GetTranslation();

    G4StepPoint* poststep = step->GetPostStepPoint();
    G4ThreeVector postPosition = poststep->GetPosition();
    G4double postKE = poststep->GetKineticEnergy();

    // Get the process name of the vertex of that particle
    if (track->GetCreatorProcess())
        creator_process_name = track->GetCreatorProcess()->GetProcessName();

    G4cout << "----------" << G4endl;
    G4cout << "Particle : " << particle_name << G4endl;
    G4cout << "stepID : " << stepID << G4endl;
    G4cout << "trackID : " << trackID << G4endl;
    G4cout << "eventID : " << eventID << G4endl;
    G4cout << "Creator Process : " << creator_process_name << G4endl;
    G4cout << "Detector name :" << physVol_name << G4endl;/*
    G4cout << "Detector position is:" << postDetectorPosition/cm << " cm" << G4endl;*/
    G4cout << "Position : " << postPosition/mm << "mm" << G4endl;
    G4cout << "Kinetic Energy is:" << postKE/MeV << " MeV" << G4endl;
    G4cout << "----------" << G4endl;
}

void Detect_edep::ClearVectorsCounts()
{
    //G4cout << "Clearing photon counts: opticalPhotonCounts size=" << opticalPhotonCounts.size() 
    //       << ", exitData size=" << exitData.size() << G4endl;
    CurrentData.clear();
    edep_per_detector.clear();
    first_time_per_detector.clear(); // New: Clear the time map
}