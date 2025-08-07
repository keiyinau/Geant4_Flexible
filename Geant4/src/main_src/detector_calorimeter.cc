#include "detector_calorimeter.hh"

Calorimeter::Calorimeter(G4String name) : G4VSensitiveDetector(name), fHitsCollectionID(-1)
{
    ClearVectorsCounts(); // Initialize the vectors to store accumulated data
	collectionName.insert("Calorimeter");


	double signalLength=500; //ns
	double SampleTime=1; //ns
	double DarkCountRate=1.7*1000*1000; //Hz
	double RiseTime=10; //ns
	double FallTimeFast=55; //ns
	double RecoveryTime=55; //ns
	double Dcr=1.7*1000*1000; //Hz
	double Xt=0.23; //ns
	double Ap=0.01; //ns	
	double pitch=40; //um
	double nCells=8334; //total number of cells
	double size=4.5; //mm
	// Reading PDE
	std::vector<double> wlen;
	std::vector<double> pde;
	std::ifstream datafile;
	datafile.open("SiPMPDE.txt");
	while(1){
		double wavelength_, pde_;
		datafile >> wavelength_ >> pde_;
		if(datafile.eof()) break; // End of file reached
		// Process the wavelength and pde values as needed
		std::cout << "Wavelength: " << wavelength_ << ", PDE: " << pde_ << std::endl;
		wlen.push_back(wavelength_);
		pde.push_back(pde_/100.0); // Convert percentage to fraction
	}
	datafile.close();
	// Electronic parameters
	double gatewidth=40; //ns
	double threshold=-0.5; //mV

	//Load the SiPM properties
	sipm::SiPMProperties myProperties = sipm::SiPMProperties();
	myProperties.setDcr(DarkCountRate);
	myProperties.setFallTimeFast(FallTimeFast);
	myProperties.setProperty("nCells",nCells);
	myProperties.setProperty("Xt",Xt);
	myProperties.setProperty("Ap",Ap);
	myProperties.setProperty("Pitch", pitch);
	myProperties.setProperty("recoveryTime", RecoveryTime);
	myProperties.setProperty("signalLength",signalLength);
	myProperties.setProperty("sampling",SampleTime);
	myProperties.setProperty("size",size);
	myProperties.setPdeType(sipm::SiPMProperties::PdeType::kSpectrumPde);
	myProperties.setPdeSpectrum(wlen,pde);
	std::cout<<"Properties:"<<myProperties<<std::endl;
	sipm::SiPMSensor mySensor(myProperties);
	std::cout<<"My Sensor:"<<mySensor<<"\n";

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
	
	mySensor.resetState();
	mySensor.addPhotons(photonTimes, photonWavelengths);
	mySensor.runEvent();

	//std::cout<<"Hits:"<<mySensor.hits().size()<<"\n";
	//std::cout<<"Debug info:"<<mySensor.debug()<<"\n";
//
	//std::cout<<"Photon count from debug info:"<<mySensor.debug().nPhotons<<"\n";
	//std::cout<<"Photonelectron count from debug info:"<<mySensor.debug().nPhotoelectrons<<"\n";
//
	//std::cout<<"Signal:"<<mySensor.signal()<<"\n";
	//SiPMAnalogSignal mySignal = mySensor.signal();
	//std::vector<float> waveform = mySignal.waveform();
	//for (float& val : waveform) {
    //        val *= -2.54;
    //    }
	SaveToRoot();
	ClearVectorsCounts(); // Clear the accumulated counts at the end of each event
}

G4bool Calorimeter::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
	
	G4Track* track = aStep->GetTrack();
	G4String particleName = track->GetParticleDefinition()->GetParticleName();
	if (particleName == "opticalphoton")
    {
		track->SetTrackStatus(fStopAndKill); // Stop the optical photon track
		SaveToStepData(aStep, ROhist, track); // Save step data for optical photons;
        return true;
    }
	return 0;
}

// This function store information to a Ntuple then it can be saved in run.cc
void Calorimeter::SaveToStepData(G4Step* aStep, G4TouchableHistory* ROhist, G4Track* track){

    G4AnalysisManager *man = G4AnalysisManager::Instance();
	G4String detector_Name = track->GetTouchable()->GetVolume()->GetName();
	G4StepPoint *preStepPoint=aStep->GetPreStepPoint();
	G4double time=preStepPoint->GetGlobalTime();
	G4ThreeVector momPhoton = preStepPoint->GetMomentum();
	G4double wlen= (1.239841939*eV/momPhoton.mag())*1E+3
	StepData data;
	data.detectorName = detector_Name;
	data.wavelength = (double)wlen;
	data.Hittime = (double)time;
	photonTimes.push_back(data.Hittime); // Store the data for this step
	photonWavelengths.push_back(data.wavelength); // Store the wavelength for this step
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
		analysisManager->FillNtupleDColumn(1,4, data.Hittime);
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
	photonTimes.clear(); // Clear the vector that stores photon times
	photonWavelengths.clear(); // Clear the vector that stores photon wavelengths
	// Clear the vector that stores the current data
	CurrentData.clear();

}