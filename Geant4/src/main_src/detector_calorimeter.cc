#include "detector_calorimeter.hh"

Calorimeter::Calorimeter(G4String name) : G4VSensitiveDetector(name), fHitsCollectionID(-1)
{
    ClearVectorsCounts(); // Initialize the vectors to store accumulated data
	collectionName.insert("Calorimeter");
}

Calorimeter::~Calorimeter()
{}

void Calorimeter::Initialize(G4HCofThisEvent* hce)
{
    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    //G4cout << "MySensitiveDetector::Initialize called for Event=" << eventID << G4endl;
    if (fHitsCollectionID < 0) {
        fHitsCollectionID = GetCollectionID(0);
    }
    G4VHitsCollection* hc = new G4VHitsCollection(SensitiveDetectorName, collectionName[0]);
    hce->AddHitsCollection(fHitsCollectionID, hc);
}

void Calorimeter::EndOfEvent(G4HCofThisEvent*){
	SaveToRoot();
	ClearVectorsCounts(); // Clear the accumulated counts at the end of each event
}

G4bool Calorimeter::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
	
	G4Track* track = aStep->GetTrack();
	G4String particleName = track->GetParticleDefinition()->GetParticleName();
	if (particleName == "opticalphoton")
    {
		SaveToStepData(aStep, ROhist, track); // Save step data for optical photons;
        return true;
    }
	return 0;
}

// This function store information to a Ntuple then it can be saved in run.cc
void Calorimeter::SaveToStepData(G4Step* aStep, G4TouchableHistory* ROhist, G4Track* track){

    G4AnalysisManager *man = G4AnalysisManager::Instance();
	G4String particle_name = track->GetDefinition()->GetParticleName();								//Get the particle name
	G4String detector_Name = track->GetTouchable()->GetVolume()->GetName();
	G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
	G4int trackID = track->GetTrackID();
	G4int parentID=aStep->GetTrack()->GetParentID();
	G4double tdep= aStep->GetDeltaTime(); // Time difference for this step
	G4StepPoint *preStepPoint=aStep->GetPreStepPoint();
	G4double time=preStepPoint->GetGlobalTime();
	scintillatorCount[parentID]+=1;
    HitTime[parentID]=time;
	StepData data;
    data.eventID = evt;
    data.trackID = trackID;
	data.detectorName = detector_Name;
	data.scintillatorCount = scintillatorCount[trackID];
	data.Hittime = HitTime[trackID];
	CurrentData.push_back(data); // Store the data for this step
	track->SetTrackStatus(fStopAndKill);

	
}
void Calorimeter::SaveToRoot(){
	std::map<G4int, StepData> filteredExitData;
    for (const auto& data : CurrentData) {
        filteredExitData[data.trackID] = data; // Overwrite with the last entry
    }
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	for (const auto&pair : filteredExitData) {
		auto& data=pair.second;
		if(data.scintillatorCount == 0) continue; // Skip if no scintillators were hit
		analysisManager->FillNtupleIColumn(1,0, data.eventID);
		analysisManager->FillNtupleIColumn(1,1, data.trackID);
		analysisManager->FillNtupleSColumn(1,2, data.detectorName);
		analysisManager->FillNtupleDColumn(1,3, data.scintillatorCount);
		analysisManager->FillNtupleDColumn(1,4, data.Hittime/ns);
		analysisManager->AddNtupleRow(1);
	}
}

// Output Information just touch the detector
void Calorimeter::ReadOut(G4Step* step, G4Track* track) {

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

void Calorimeter::ClearVectorsCounts()
{
    //G4cout << "Clearing photon counts: opticalPhotonCounts size=" << opticalPhotonCounts.size() 
    //       << ", exitData size=" << exitData.size() << G4endl;
	scintillatorCount.clear();
	HitTime.clear();
	// Clear the vector that stores the current data
	CurrentData.clear();

}