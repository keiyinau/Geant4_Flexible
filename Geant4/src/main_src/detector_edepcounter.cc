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
    edep=0.;
	G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    //G4cout << "MySensitiveDetector::Initialize called for Event=" << eventID << G4endl;
    if (fHitsCollectionID < 0) {
        fHitsCollectionID = GetCollectionID(0);
    }
    G4VHitsCollection* hc = new G4VHitsCollection(SensitiveDetectorName, collectionName[0]);
    hce->AddHitsCollection(fHitsCollectionID, hc);
}

void Detect_edep::EndOfEvent(G4HCofThisEvent*){
	SaveToRoot();
	ClearVectorsCounts(); // Clear the accumulated counts at the end of each event
}

G4bool Detect_edep::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
	
	G4Track* track = aStep->GetTrack();
	SaveToStepData(aStep,ROhist,track);

	return 0;
}

// This function store information to a Ntuple then it can be saved in run.cc
void Detect_edep::SaveToStepData(G4Step* aStep, G4TouchableHistory* ROhist, G4Track* track){
	edep += aStep->GetTotalEnergyDeposit();
	
    G4AnalysisManager *man = G4AnalysisManager::Instance();
	G4String detector_Name = track->GetTouchable()->GetVolume()->GetName();
	G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
	
	data.eventID = evt;
	data.detectorName = detector_Name;
	data.edep_accumulated = edep;
}
void Detect_edep::SaveToRoot(){
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->FillNtupleIColumn(1, 0, data.eventID); // eventID
	analysisManager->FillNtupleSColumn(1, 1, data.detectorName); // trackID
	analysisManager->FillNtupleDColumn(1, 2, data.edep_accumulated/MeV); // stepID
	analysisManager->AddNtupleRow(1);
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

}