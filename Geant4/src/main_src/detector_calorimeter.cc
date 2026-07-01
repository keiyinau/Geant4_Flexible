#include "detector_calorimeter.hh"
#include "event.hh"
Calorimeter::Calorimeter(G4String name) : G4VSensitiveDetector(name), fHitsCollectionID(-1)
{
    ClearVectorsCounts(); // Initialize the vectors to store accumulated data
	collectionName.insert("Calorimeter");
	isGraph=false;
	isDCR=true;
	isXT=true;
	isAP=true;
	signalLength=1000; //ns
	SampleTime=1; //ns
	DarkCountRate=2.04*1e6; //Hz
	RiseTime=10; //ns
	FallTimeFast=200; //ns
	RecoveryTime=55; //ns
	Dcr=2.04*1e6; //Hz
	Xt=0.25; //ns
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
	gatewidth=250; //ns
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

    // Get primary decay time from event action
    const MyEventAction* eventAction = static_cast<const MyEventAction*>(G4RunManager::GetRunManager()->GetUserEventAction());
    G4double primaryTime = eventAction->GetPrimaryDecayTime();

    for (const auto& pair : photonTimes_per_detector) {
        G4String det_name = pair.first;
        const auto& times = pair.second;
        const auto& wlens = photonWavelengths_per_detector[det_name];

        if (times.empty()) continue;

        std::vector<double> shiftedTimes;
        shiftedTimes.reserve(times.size());
        for (double t : times) {
            double t_rel = t - primaryTime;
            shiftedTimes.push_back(t_rel);
        }

        mySensor.resetState();
        mySensor.addPhotons(shiftedTimes, wlens);
        mySensor.runEvent();

        const auto& debug = mySensor.debug();
        const auto& signal = mySensor.signal();

        // === Find first REAL photon time (from CsI) ===
        G4double firstPhotonTime = -1.0;
        if (!shiftedTimes.empty()) {
            firstPhotonTime = *std::min_element(shiftedTimes.begin(), shiftedTimes.end());
        }

        // === Only process if at least one real photon ===
        if (firstPhotonTime >= 0) {
            G4double gateStart = firstPhotonTime;
            G4double gateEnd   = gateStart + gatewidth;

            LoadData data;
            data.eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
            data.detectorName = det_name;
            data.RealPhotonCount = debug.nPhotons;
            data.PEsCount = debug.nPhotoelectrons;
            data.NoisePEsCount = debug.nDcr + debug.nXt + debug.nAp;
            data.Time_Of_Triggering = firstPhotonTime;  // Real trigger time

            // Ensure gate fits in signal length
            if (gateEnd <= signalLength) {
                G4double integral = signal.integral(gateStart, gateEnd, 14.0);  // No threshold for integration

                if (integral < 1e10) {  // Removed >0 to save even if integral==0
                    data.Area = integral;
                    CurrentData.push_back(data);
                }
            }
            else{
                data.Area = -100.0; // Indicate invalid area due to gate exceeding signal length
                CurrentData.push_back(data);
            }
        }

        // Optional: Plot waveform
        if (isGraph) {
            PlotWaveform(signal, det_name);
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
	photonTimes_per_detector[detector_Name].push_back(time);
	photonWavelengths_per_detector[detector_Name].push_back(wlen);
}
void Calorimeter::SaveToRoot(){
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	const G4double SiPM_gain=1e6;
	const G4double preamp_gain=20.0;
	const G4double electron_charge=1.6217662e-19;
	const G4double ADC_gain=100*1e-15;
	const G4double Total_Gain=ADC_gain/(preamp_gain*SiPM_gain*electron_charge);
	const G4double unit_Area=163.73262320291757;
	for(const auto&data:CurrentData){
		//calibrated_pe=data.PEsCount;
		//calibrated_area=data.Area;
		analysisManager->FillNtupleIColumn(0,0, data.eventID);
		analysisManager->FillNtupleSColumn(0,1, data.detectorName);
		analysisManager->FillNtupleDColumn(0,2, data.Area);
		analysisManager->FillNtupleIColumn(0,3, data.RealPhotonCount);
		analysisManager->FillNtupleIColumn(0,4, data.PEsCount);
		analysisManager->FillNtupleIColumn(0,5, data.NoisePEsCount);
		analysisManager->FillNtupleDColumn(0,6, data.Time_Of_Triggering);
		analysisManager->FillNtupleDColumn(0,7, data.Area/unit_Area);

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
