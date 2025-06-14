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
	if(particleName == "gamma" || particleName == "e-") {
		SaveToStepData(aStep,ROhist,track);
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
	AccumatedDistance_count[trackID]+=ldep;
    AccumulatedTime_count[trackID]+=tdep;
    AccumulatedEnergy_count[trackID]+=edep;
	StepData data;
    data.eventID = evt;
    data.trackID = trackID;
	data.stepID = stepID;
	data.particleName = particle_name;
	data.ProcessName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName(); // Get the process name of the step
	data.kineticEnergy = ekin;
	data.AccumatedDistance = AccumatedDistance_count[trackID];
	data.AccumulatedTime = AccumulatedTime_count[trackID];
	data.AccumulatedEnergy = AccumulatedEnergy_count[trackID];
	CurrentData.push_back(data); // Store the data for this step

	
}
void Tracker::SaveToRoot(){
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	for(const auto&data:CurrentData){
		// Fill the ntuple with the data
		analysisManager->FillNtupleIColumn(0,0, data.eventID);
		analysisManager->FillNtupleIColumn(0,1, data.trackID);
		analysisManager->FillNtupleIColumn(0,2, data.stepID);
		analysisManager->FillNtupleSColumn(0,3, data.particleName);
		analysisManager->FillNtupleSColumn(0,4, data.ProcessName); // Fill the process name
		analysisManager->FillNtupleDColumn(0,5, data.kineticEnergy/MeV);
		analysisManager->FillNtupleDColumn(0,6, data.AccumatedDistance/mm); // Fill the accumulated distance
		analysisManager->FillNtupleDColumn(0,7, data.AccumulatedTime/ns); // Fill the accumulated time
		analysisManager->FillNtupleDColumn(0,8, data.AccumulatedEnergy/MeV); // Fill the accumulated energy
		analysisManager->AddNtupleRow(0);
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