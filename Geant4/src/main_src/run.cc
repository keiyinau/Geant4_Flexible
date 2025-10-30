#include "run.hh"
#include "construction.hh"

MyRunAction::MyRunAction(){
	G4AnalysisManager *man = G4AnalysisManager::Instance();

    // Here to select which type of Data File to be created. Options: Normal, Vertex, Step
    G4int select = 0;
    if (select == 0)
        CreateDataFile_SensitiveDetector(man);
}

MyRunAction::~MyRunAction()
{}

void MyRunAction::BeginOfRunAction(const G4Run* run){
	G4AnalysisManager* man = G4AnalysisManager::Instance();
	man->SetNtupleMerging(true);
	G4int runID = run->GetRunID();
    G4String file_name = MyDetectorConstruction::file_name;

	std::stringstream strRunID;
	strRunID << runID;

    if (file_name == "")
	    man->OpenFile("output_"+strRunID.str()+".root");
    else
        man->OpenFile(file_name+".root");
}

void MyRunAction::EndOfRunAction(const G4Run* run){
	G4AnalysisManager* man = G4AnalysisManager::Instance();
	man->Write();
	man->CloseFile();
}

void MyRunAction::CreateDataFile_SensitiveDetector(G4AnalysisManager* man)
{
    man->CreateNtuple("SiPM", "SiPM Data"); // In practice, only X,Y,Z,Edep,time,track is available
    man->CreateNtupleIColumn("EventID");
    man->CreateNtupleDColumn("Area"); // Accumulated distance in mm
    man->CreateNtupleIColumn("#RealPhoton"); // Accumulated
    man->CreateNtupleIColumn("#PE"); // Accumulated
    man->CreateNtupleIColumn("#NoisePE"); // Position X in
    man->CreateNtupleDColumn("Time_Of_Triggering_ns"); // Position Y in
    //man->CreateNtupleSColumn("SiPMName");
    man->FinishNtuple(0);

    man->CreateNtuple("Ideal", "Ideal data"); 
    man->CreateNtupleIColumn("EventID");
    man->CreateNtupleIColumn("TrackID");
    man->CreateNtupleIColumn("StepID");
    man->CreateNtupleIColumn("ParentID");
    man->CreateNtupleSColumn("DetectorName");
    man->CreateNtupleSColumn("ParticleName");
    man->CreateNtupleSColumn("CreatorProcessName");
    man->CreateNtupleSColumn("ProcessName");
    man->CreateNtupleDColumn("KineticEnergy_MeV"); // in Me
    man->CreateNtupleDColumn("x_distance_mm");
    man->CreateNtupleDColumn("y_distance_mm");
    man->CreateNtupleDColumn("z_distance_mm");
    man->CreateNtupleDColumn("x_momentum_MeV");
    man->CreateNtupleDColumn("y_momentum_MeV");
    man->CreateNtupleDColumn("z_momentum_MeV");
    man->FinishNtuple(1);


    man->CreateNtuple("Tracker", "Tracker data"); 
    man->CreateNtupleIColumn("EventID");
    man->CreateNtupleIColumn("TrackID");
    man->CreateNtupleDColumn("AccumatedDistance_mm");
    man->CreateNtupleDColumn("AccumulatedTime_ns");
    man->CreateNtupleDColumn("AccumulatedEnergy_MeV");
    man->CreateNtupleSColumn("DetectorName");
    man->CreateNtupleSColumn("ParticleName");
    man->CreateNtupleSColumn("CreatorProcessName");
    man->FinishNtuple(2);

}
