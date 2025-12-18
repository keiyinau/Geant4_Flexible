#include "detector_calorimeter.hh"

Calorimeter::Calorimeter(G4String name) : G4VSensitiveDetector(name), fHitsCollectionID(-1)
{
    ClearVectorsCounts(); // Initialize the vectors to store accumulated data
	collectionName.insert("Calorimeter");
	isGraph=false;
	isDCR=false;
	isXT=true;
	isAP=true;
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
		//std::cout << "Wavelength: " << wavelength_ << ", PDE: " << pde_ << std::endl;
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
	if(isXT==false){
		myProperties.setXtOff();
	}
	if(isDCR==false){
		myProperties.setDcrOff();
	}
	if(isAP==false){
		myProperties.setApOff();
	}
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
    if (fHitsCollectionID < 0) {
        fHitsCollectionID = GetCollectionID(0);
    }
    G4VHitsCollection* hc = new G4VHitsCollection(SensitiveDetectorName, collectionName[0]);
    hce->AddHitsCollection(fHitsCollectionID, hc);
	
	
}

void Calorimeter::EndOfEvent(G4HCofThisEvent*)
{
	if (photonTimes_per_detector.empty()) {
        ClearVectorsCounts();
        return;
    }

    // Compute global minimum time across all detectors with photons
    G4double global_min_time = DBL_MAX;
    for (const auto& pair : photonTimes_per_detector) {
        const auto& times = pair.second;
        if (!times.empty()) {
            double min_t = *std::min_element(times.begin(), times.end());
            if (min_t < global_min_time) {
                global_min_time = min_t;
            }
        }
    }
    if (global_min_time == DBL_MAX) {
        global_min_time = 0.;
    }

    // Process each detector separately
    for (const auto& pair : photonTimes_per_detector) {
        G4String det_name = pair.first;
        const auto& times = pair.second;
        const auto& wlens = photonWavelengths_per_detector[det_name];

        if (times.empty()) continue;

        // Compute per-detector min time for relative triggering
        double det_min_time = *std::min_element(times.begin(), times.end());

        // Shift times relative to det_min_time for sensor processing
        std::vector<double> shiftedTimes;
        shiftedTimes.reserve(times.size());
        for (double t : times) {
            double t_rel = t - det_min_time;
            shiftedTimes.push_back(t_rel);
        }

        mySensor.resetState();
        mySensor.addPhotons(shiftedTimes, wlens);
        mySensor.runEvent();

        const auto& debug = mySensor.debug();
        const auto& signal = mySensor.signal();

        // Relative triggering time: det_min_time - global_min_time
        G4double rel_trigger_time = (det_min_time - global_min_time) / ns;  // in ns

        // Gate starts at 0 in relative time (first photon in this det)
        G4double gateStart = 0.0;
        G4double gateEnd = gateStart + gatewidth;

        if (gateEnd <= signalLength) {
            G4double integral = signal.integral(gateStart, gateEnd, 0.0);

            if (integral < 1e10) {
                LoadData data;
                data.eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
                data.detectorName = det_name;  // New: store detector name
                data.Area = integral;
                data.RealPhotonCount = debug.nPhotons;
                data.PEsCount = debug.nPhotoelectrons;
                data.NoisePEsCount = debug.nDcr + debug.nXt + debug.nAp;
                data.Time_Of_Triggering = rel_trigger_time;  // Relative time in ns

                CurrentData.push_back(data);
            }
        } else {
            LoadData data;
            data.eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
            data.detectorName = det_name;  // New: store detector name
            data.Area = -100.0;
            data.RealPhotonCount = debug.nPhotons;
            data.PEsCount = debug.nPhotoelectrons;
            data.NoisePEsCount = debug.nDcr + debug.nXt + debug.nAp;
            data.Time_Of_Triggering = rel_trigger_time;
            CurrentData.push_back(data);
        }

        // Optional: Plot waveform per detector
        if (isGraph) {
            PlotWaveform(signal, det_name);  // Modify PlotWaveform to take det_name for unique filename
        }
    }

    SaveToRoot();
    ClearVectorsCounts();
}

void Calorimeter::PlotWaveform(const sipm::SiPMAnalogSignal& signal, const G4String& det_name)
{
    std::vector<float> waveform = signal.waveform();
    for (float& val : waveform) val *= gain;

    size_t nPoints = waveform.size();
    TGraph* graph = new TGraph(nPoints);
    for (size_t i = 0; i < nPoints; ++i) {
        graph->SetPoint(i, i * SampleTime, waveform[i]);
    }

    TCanvas* c = new TCanvas("c", "SiPM Signal", 800, 600);
    graph->SetTitle(("SiPM Waveform - " + det_name + ";Time (ns);Amplitude (mV)").c_str());
    graph->Draw("AL");
    c->SaveAs(("waveform_" + det_name + ".png").c_str());

    delete graph;
    delete c;
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
	//StepData data;  // No longer needed
	//data.detector_Name = detector_Name;
	//data.wavelength = (double)wlen;
	//data.Hittime = (double)time;
	photonTimes_per_detector[detector_Name].push_back(time); // Store per detector
	photonWavelengths_per_detector[detector_Name].push_back(wlen); // Store per detector
	//G4cout<<"Photon Hit at"<<data.Hittime/ns<<std::endl;
	//std::cout<<wlen<<std::endl;
}
void Calorimeter::SaveToRoot(){
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	for(const auto&data:CurrentData){
		analysisManager->FillNtupleIColumn(0,0, data.eventID);
		analysisManager->FillNtupleSColumn(0,1, data.detectorName); // Uncommented and using detectorName
		analysisManager->FillNtupleDColumn(0,2, data.Area);
		analysisManager->FillNtupleIColumn(0,3, data.RealPhotonCount);
		analysisManager->FillNtupleIColumn(0,4, data.PEsCount);
		analysisManager->FillNtupleIColumn(0,5, data.NoisePEsCount);
		analysisManager->FillNtupleDColumn(0,6, data.Time_Of_Triggering);
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
	photonTimes_per_detector.clear(); // Clear the map that stores photon times per detector
	photonWavelengths_per_detector.clear(); // Clear the map that stores photon wavelengths per detector
	detectorname.clear();
	// Clear the vector that stores the current data
	CurrentData.clear();

}