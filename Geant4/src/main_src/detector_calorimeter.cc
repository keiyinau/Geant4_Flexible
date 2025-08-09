#include "detector_calorimeter.hh"

Calorimeter::Calorimeter(G4String name) : G4VSensitiveDetector(name), fHitsCollectionID(-1)
{
    ClearVectorsCounts(); // Initialize the vectors to store accumulated data
	collectionName.insert("Calorimeter");
	isGraph=false;
	signalLength=500; //ns
	SampleTime=1; //ns
	DarkCountRate=1.7*1000*1000; //Hz
	RiseTime=10; //ns
	FallTimeFast=55; //ns
	RecoveryTime=55; //ns
	Dcr=1.7*1000*1000; //Hz
	Xt=0.23; //ns
	Ap=0.01; //ns	
	pitch=40; //um
	nCells=8334; //total number of cells
	size=4.5; //mm
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
	gatewidth=40; //ns
	threshold=0.5; //mV
	gain=-2.5;

	//Load the SiPM properties
	myProperties = sipm::SiPMProperties();
	myProperties.setDcr(DarkCountRate);
	myProperties.setFallTimeFast(FallTimeFast);
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
	mySensor = sipm::SiPMSensor(myProperties);
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
	//for(int i=0; i<photonTimes.size(); i++){
	//	//std::cout<<"Photon "<<i<<": Time = "<<photonTimes[i]<<" ns, Wavelength = "<<photonWavelengths[i]<<" nm"<<std::endl;
	//	std::cout<<""<<photonTimes[i]<<std::endl;
	//}
	//std::cout<<"Debug info:"<<mySensor.debug()<<"\n";
	//std::cout<<"Photon count from debug info:"<<mySensor.debug().nPhotons<<"\n";
	//std::cout<<"Photonelectron count from debug info:"<<mySensor.debug().nPhotoelectrons<<"\n";
	//std::cout<<"Number of DCR event"<<mySensor.debug().nDcr<<"\n";
	//std::cout<<"Number of Optical Cross talk event"<<mySensor.debug().nXt<<"\n";
	//std::cout<<"Number of After pulse event"<<mySensor.debug().nAp<<"\n";
//
	//std::cout<<"Signal:"<<mySensor.signal()<<"\n";
	mySignal = mySensor.signal();
	G4double timeofArrival= mySignal.toa(0,signalLength, threshold);
  	G4double integral = gain*mySignal.integral((double)timeofArrival,(double)timeofArrival+(double)gatewidth,threshold);
	LoadData data;
	data.eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
	data.Area = integral;
	data.RealPhotonCount = mySensor.debug().nPhotons;
	data.PEsCount = mySensor.debug().nPhotoelectrons;
	data.NoisePEsCount = mySensor.debug().nDcr + mySensor.debug().nXt + mySensor.debug().nAp;
	data.Time_Of_Triggering = timeofArrival;
	CurrentData.push_back(data); // Store the data for this event


	if(isGraph){
		std::vector<float> waveform = mySignal.waveform();
		for (float& val : waveform) {
				val *= gain;
		}


		// Generate time points for plotting
		size_t nPoints = static_cast<size_t>(signalLength / SampleTime);
		std::vector<double> timePoints(nPoints);
		for (size_t i = 0; i < nPoints; ++i) {
			timePoints[i] = i * SampleTime;
		}
		// Create ROOT TGraph for plotting
		TGraph* graph = new TGraph(nPoints);
		for (size_t i = 0; i < nPoints && i < waveform.size(); ++i) {
			graph->SetPoint(i, timePoints[i], waveform[i]);
		}
		graph->SetTitle("SiPM Waveform;Time (ns);Amplitude (mV)");
		graph->SetMarkerStyle(20);
		graph->SetMarkerSize(0.5);

		// Create canvas and draw the graph
		TCanvas* canvas = new TCanvas("canvas", "SiPM Waveform", 800, 600);
		graph->Draw("APL"); // Draw with axis, points, and line
		canvas->Update();

		// Save the plot to a file
		canvas->SaveAs("waveform.png");

		// If running interactively, uncomment the next line to keep the window open
		// app.Run();

		// Clean up
		delete graph;
		delete canvas;
	}
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
	G4double wlen= (1239.841939/(track->GetDynamicParticle()->GetTotalEnergy()/eV));
	StepData data;
	data.detector_Name = detector_Name;
	data.wavelength = (double)wlen;
	data.Hittime = (double)time;
	photonTimes.push_back(data.Hittime); // Store the data for this step
	photonWavelengths.push_back(data.wavelength); // Store the wavelength for this step
	//G4cout<<"Photon Hit at"<<data.Hittime/ns<<std::endl;
	//std::cout<<wlen<<std::endl;
}
void Calorimeter::SaveToRoot(){
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	for(const auto&data:CurrentData){
		analysisManager->FillNtupleIColumn(0,0, data.eventID);
		//analysisManager->FillNtupleSColumn(0,1, data.SiPMName);
		analysisManager->FillNtupleDColumn(0,1, data.Area);
		analysisManager->FillNtupleIColumn(0,2, data.RealPhotonCount);
		analysisManager->FillNtupleIColumn(0,3, data.PEsCount);
		analysisManager->FillNtupleIColumn(0,4, data.NoisePEsCount);
		analysisManager->FillNtupleDColumn(0,5, data.Time_Of_Triggering);
		// Fill the ntuple with the data
		analysisManager->AddNtupleRow(0);
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
	detectorname.clear();
	// Clear the vector that stores the current data
	CurrentData.clear();

}