#include "detector_tracker.hh"

Tracker::Tracker(G4String name) : G4VSensitiveDetector(name), fHitsCollectionID(-1)
{
    ClearVectorsCounts(); // Initialize the vectors to store accumulated data
	collectionName.insert("Tracker");
}

Tracker::~Tracker()
{}

void Tracker::Initialize(G4HCofThisEvent* hce)
{
    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    //G4cout << "MySensitiveDetector::Initialize called for Event=" << eventID << G4endl;
    if (fHitsCollectionID < 0) {
        fHitsCollectionID = GetCollectionID(0);
    }
    G4VHitsCollection* hc = new G4VHitsCollection(SensitiveDetectorName, collectionName[0]);
    hce->AddHitsCollection(fHitsCollectionID, hc);
}

void Tracker::EndOfEvent(G4HCofThisEvent*){
	SaveToRoot();
	ClearVectorsCounts(); // Clear the accumulated counts at the end of each event
}

G4bool Tracker::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
	G4Track* track = aStep->GetTrack();
	G4String particleName = track->GetParticleDefinition()->GetParticleName(); //Consider only Positron

    if (particleName == "opticalphoton") {
        return false;
    }
    if (track->GetParentID() == 0 || particleName == "gamma") {
        SaveToStepData(aStep, ROhist, track);
        return true;
    }
	return 0;
}

// This function store information to a Ntuple then it can be saved in run.cc
void Tracker::SaveToStepData(G4Step* aStep, G4TouchableHistory* ROhist, G4Track* track){

    G4AnalysisManager *man = G4AnalysisManager::Instance();
	G4String particle_name = track->GetDefinition()->GetParticleName();								//Get the particle name
	G4String detector_Name = track->GetTouchable()->GetVolume()->GetName();
	G4String creator_process_name = "NULL";
	if (track->GetCreatorProcess())
		creator_process_name = track->GetCreatorProcess()->GetProcessName();				//Get the process name of the vertex of track 
	G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
	G4int trackID = track->GetTrackID();
	G4int stepID = track->GetCurrentStepNumber();
	G4int parentID=aStep->GetTrack()->GetParentID();
	G4double edep = aStep->GetTotalEnergyDeposit(); // Total energy deposited in this step
	G4double ldep= aStep->GetStepLength(); // Length of the step
	G4double tdep= aStep->GetDeltaTime(); // Time difference for this step
	G4double ekin = track->GetKineticEnergy();	
	G4String particlename = track->GetDefinition()->GetParticleName(); // Get the particle name
	AccumatedDistance_count[trackID]+=ldep;
    AccumulatedTime_count[trackID]+=tdep;
    AccumulatedEnergy_count[trackID]+=edep;
	StepData data;
    data.eventID = evt;
    data.trackID = trackID;
	data.AccumatedDistance = AccumatedDistance_count[trackID];
	data.AccumulatedTime = AccumulatedTime_count[trackID];
	data.AccumulatedEnergy = AccumulatedEnergy_count[trackID];
	data.postPositionX = aStep->GetPostStepPoint()->GetPosition().x();
	data.postPositionY = aStep->GetPostStepPoint()->GetPosition().y();
	data.postPositionZ = aStep->GetPostStepPoint()->GetPosition().z();
	if(StepData_count[trackID].preKE<aStep->GetPreStepPoint()->GetKineticEnergy()){ // Get the kinetic energy at the pre-step point) { // If this is the first step for this track
		data.preKE = aStep->GetPreStepPoint()->GetKineticEnergy(); // Set the event ID
	} else {
		data.preKE = StepData_count[trackID].preKE; // Use the existing event ID
	}
	data.detectorName = detector_Name; // Store the detector name
	data.particleName = particle_name; // Store the particle name
	data.creatorProcessName = creator_process_name; // Store the creator process name
	StepData_count[trackID] = data; // Store the data in the map
}
void Tracker::SaveToRoot(){
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	for (const auto& pair : StepData_count) {
		const auto& data = pair.second; // Access the StepData object
		//G4cout<<data.AccumatedDistance / mm;
		analysisManager->FillNtupleIColumn(2, 0, data.eventID);
		analysisManager->FillNtupleIColumn(2, 1, data.trackID);
		analysisManager->FillNtupleDColumn(2, 2, data.AccumatedDistance / mm);
		analysisManager->FillNtupleDColumn(2, 3, data.AccumulatedTime / ns);
		analysisManager->FillNtupleDColumn(2, 4, data.AccumulatedEnergy / MeV);
		analysisManager->FillNtupleSColumn(2, 5, data.detectorName);
		analysisManager->FillNtupleSColumn(2, 6, data.particleName);
		analysisManager->FillNtupleSColumn(2, 7, data.creatorProcessName);
		analysisManager->AddNtupleRow(2);
	}
}

// Output Information just touch the detector
void Tracker::ReadOut(G4Step* step, G4Track* track) {

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

void Tracker::ClearVectorsCounts()
{
    //G4cout << "Clearing photon counts: opticalPhotonCounts size=" << opticalPhotonCounts.size() 
    //       << ", exitData size=" << exitData.size() << G4endl;
    AccumatedDistance_count.clear();
    AccumulatedTime_count.clear();
    AccumulatedEnergy_count.clear();
	CurrentData.clear();

}