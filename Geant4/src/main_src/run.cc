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
    man->CreateNtuple("Tracker", "Tracker Data"); // In practice, only X,Y,Z,Edep,time,track is available
    man->CreateNtupleIColumn("EventID");
    //man->CreateNtupleIColumn("TrackID");
    man->CreateNtupleDColumn("AccumulatedDistance_mm"); // Accumulated distance in mm
    //man->CreateNtupleDColumn("AccumulatedTime_ns"); // Accumulated
    //man->CreateNtupleDColumn("AccumulatedEnergy_MeV"); // Accumulated
    //man->CreateNtupleDColumn("PositionX_mm"); // Position X in
    //man->CreateNtupleDColumn("PositionY_mm"); // Position Y in
    man->CreateNtupleDColumn("PositionZ_mm"); // Position Z in
    man->CreateNtupleDColumn("KineticEnergy_MeV"); // Kinetic Energy
    //man->CreateNtupleSColumn("DetectorName"); // Detector Name
    //man->CreateNtupleSColumn("ParticleName"); // Particle Name
    //man->CreateNtupleSColumn("CreatorProcess"); // Creator Process
    man->FinishNtuple(0);
//
    //man->CreateNtuple("Reference", "Reference"); // Only record the enter, leaving datas
    //man->CreateNtupleIColumn("EventID");
    //man->CreateNtupleIColumn("TrackID");
    //man->CreateNtupleIColumn("StepID");
    //man->CreateNtupleIColumn("ParentID");
    //man->CreateNtupleSColumn("DetectorName");
    //man->CreateNtupleSColumn("ParticleName");
    //man->CreateNtupleSColumn("CreatorProcess");
    //man->CreateNtupleSColumn("ProcessName");
    //man->CreateNtupleDColumn("KineticEnergy_MeV");
    //man->CreateNtupleDColumn("PositionX_mm"); // Position X in
    //man->CreateNtupleDColumn("PositionY_mm"); // Position Y in
    //man->CreateNtupleDColumn("PositionZ_mm"); // Position Z in
    //man->FinishNtuple(1);


}
