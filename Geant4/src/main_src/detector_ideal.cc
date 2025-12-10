#include "detector_ideal.hh"

Detect_reference::Detect_reference(G4String name) : G4VSensitiveDetector(name), fHitsCollectionID(-1)
{
    ClearVectorsCounts(); // Initialize the vectors to store accumulated data
	collectionName.insert("Reference");
}

Detect_reference::~Detect_reference()
{}

void Detect_reference::Initialize(G4HCofThisEvent* hce)
{
    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    //G4cout << "MySensitiveDetector::Initialize called for Event=" << eventID << G4endl;
    if (fHitsCollectionID < 0) {
        fHitsCollectionID = GetCollectionID(0);
    }
    G4VHitsCollection* hc = new G4VHitsCollection(SensitiveDetectorName, collectionName[0]);
    hce->AddHitsCollection(fHitsCollectionID, hc);
}

void Detect_reference::EndOfEvent(G4HCofThisEvent*){
	SaveToRoot();
	ClearVectorsCounts(); // Clear the accumulated counts at the end of each event
}

G4bool Detect_reference::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
	
	G4Track* track = aStep->GetTrack();
		//track->SetTrackStatus(fStopAndKill);
	G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
	G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
	G4String particleName = track->GetParticleDefinition()->GetParticleName();
	if(particleName != "opticalphoton") {
		// Check if particle enters the detector (first entry)
		if (preStepPoint->GetStepStatus() == fGeomBoundary|| track->GetCurrentStepNumber() == 1){
			SaveToStepData(aStep,ROhist,track);
			//ReadOut(aStep, track);	// Output Information just touch the detector
			//track->SetTrackStatus(fStopAndKill);
		}

		// Check if particle leaves the detector or dies inside
		if ((postStepPoint->GetStepStatus() == fGeomBoundary || track->GetTrackStatus() == fStopAndKill)){
			// Record leave/die data
			SaveToStepData(aStep,ROhist,track);
			//ReadOut(aStep, track);	// Output Information just touch the detector
			//track->SetTrackStatus(fStopAndKill);
		}
	}
	return 0;
}

// This function store information to a Ntuple then it can be saved in run.cc
void Detect_reference::SaveToStepData(G4Step* aStep, G4TouchableHistory* ROhist, G4Track* track){
    G4AnalysisManager *man = G4AnalysisManager::Instance();
	G4String particle_name = track->GetDefinition()->GetParticleName();								//Get the particle name
	G4String detector_Name = track->GetTouchable()->GetVolume()->GetName();
	G4String creator_process_name = "NULL";
	if (track->GetCreatorProcess())
		creator_process_name = track->GetCreatorProcess()->GetProcessName();				//Get the process name of the vertex of track 
	G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
	G4int trackID = track->GetTrackID();
	G4int parentID=aStep->GetTrack()->GetParentID();
	G4double ekin = track->GetKineticEnergy();
	G4double x_pos = aStep->GetPreStepPoint()->GetPosition().x();
	G4double y_pos = aStep->GetPreStepPoint()->GetPosition().y();
	G4double z_pos = aStep->GetPreStepPoint()->GetPosition().z();	
	StepData data;
	data.eventID = evt;
	data.trackID = trackID;
	data.parentID = parentID; // Store the parent ID
	data.detectorName = detector_Name;
	data.particleName = particle_name;
	data.creatorProcessName = creator_process_name; // Process that created the track
	data.ProcessName = track->GetCreatorProcess() ? track->GetCreatorProcess()->GetProcessName() : "NULL"; // Process name
	data.kineticEnergy = ekin;
	data.x_distance = x_pos; // Store the x position
	data.y_distance = y_pos; // Store the y position
	data.z_distance = z_pos; // Store the z position
	CurrentData.push_back(data); // Store the data for this step
}
void Detect_reference::SaveToRoot(){
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	for(const auto&data:CurrentData){
		// Fill the ntuple with the data
		analysisManager->FillNtupleIColumn(3, 0, data.eventID); // eventID
		analysisManager->FillNtupleIColumn(3, 1, data.trackID); // trackID
		analysisManager->FillNtupleIColumn(3, 2, data.parentID); // parentID
		analysisManager->FillNtupleSColumn(3, 3, data.detectorName); // detectorName
		analysisManager->FillNtupleSColumn(3, 4, data.particleName); // particleName
		analysisManager->FillNtupleSColumn(3, 5, data.creatorProcessName); // creatorProcessName
		analysisManager->FillNtupleSColumn(3, 6, data.ProcessName); // ProcessName
		analysisManager->FillNtupleDColumn(3, 7, data.kineticEnergy/MeV); // kineticEnergy
		analysisManager->FillNtupleDColumn(3, 8, data.x_distance/mm); // x_distance
		analysisManager->FillNtupleDColumn(3, 9, data.y_distance/mm); // y_distance
		analysisManager->FillNtupleDColumn(3, 10, data.z_distance/mm); // z_distance
		analysisManager->AddNtupleRow(3);
	}
}

// Output Information just touch the detector
void Detect_reference::ReadOut(G4Step* step, G4Track* track) {

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

void Detect_reference::ClearVectorsCounts()
{
    //G4cout << "Clearing photon counts: opticalPhotonCounts size=" << opticalPhotonCounts.size() 
    //       << ", exitData size=" << exitData.size() << G4endl;
	CurrentData.clear();

}