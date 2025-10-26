#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include <cmath>

#include "G4VUserDetectorConstruction.hh"

#include "globals.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4GenericMessenger.hh"

#include "G4UnionSolid.hh"

#include "G4SDManager.hh"

#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "detector.hh"
#include "detector_tracker.hh"
#include "detector_ideal.hh"
#include "detector_calorimeter.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
	MyDetectorConstruction();
	~MyDetectorConstruction();

	virtual G4VPhysicalVolume* Construct();
	virtual void ConstructSDandField();
	void DefineMaterials();
	void DefineMessenger();
	bool readAndProcessData(const std::string& filename, 
				   std::vector<double>& emission_Energy, 
				   std::vector<double>& emission_fractions);
	bool readAndProcessData_Energy(const std::string& filename, 
					   std::vector<double>& emission_Energy, 
					   std::vector<double>& emission_fractions);
	bool readAndProcessData_txt(const std::string& filename, 
				   std::vector<double>& emission_Energy, 
				   std::vector<double>& emission_fractions);
	bool readAndProcessData_Energy_txt(const std::string& filename, 
			   std::vector<double>& emission_Energy, 
			   std::vector<double>& emission_fractions);
	bool readAndProcessData_Energy_cm_txt(const std::string& filename, 
		   std::vector<double>& emission_Energy, 
		   std::vector<double>& emission_fractions);
	
	// Ideal Detector
	void ConstructShell_Detector();
	// End Ideal Detector
	// tpc
	void ConstructTPC();
	// End tpc
	// Calorimeter
	void ConstructCalorimeter();
	void ConstructCalorimeter_unit(G4ThreeVector translation, G4double angle,G4String name);
	// End Calorimeter
	//Construct source
	void ConstructSource();
	void ConstructLiquidScintillator();
	// End Construct source

	static G4String file_name;

private:
	G4GenericMessenger* fMessenger;
	// World
	G4Material *Air, *Vacuum, *matWorld;
	G4LogicalVolume *logicWorld;
	G4VPhysicalVolume *physWorld;
	// End World

	// Ideal Detector
	G4bool isDetector_Shell;
	G4LogicalVolume *logicDetector_Shell;
	G4VPhysicalVolume *physDetector_Shell;
	// End Ideal Detector

	// Ideal Detector
	G4bool isTPC;
	G4LogicalVolume *logicTPC;
	G4VPhysicalVolume *physTPC;
	// End Ideal Detector
	// Ideal Detector
	G4bool isCalorimeter;
	G4LogicalVolume *logicCalorimeter;
	G4VPhysicalVolume *physCalorimeter;
	G4OpticalSurface *surfCsI_SiPM, *surfCsI_Teflon;
	std::vector<G4LogicalVolume*> logicScintillators,logicSiPM,logicTapflon;
	// End Ideal Detector
	

	// Radioactive Source (Positron Source)
	G4bool isSource;
	G4double ring_radius, ring_height_half, disk_radius, disk_height_half, bare_source_radius, bare_source_height_half;
	G4LogicalVolume *logicRing, *logicDisk, *logicBareSource;
	G4VPhysicalVolume *physRing, *physDisk, *physBareSource;
	G4bool isBareSource_Dt, isBS_Disk_Dt, isBSD_Ring_Dt;						//BS = Bare Source, D = Disk, Dt = Detector
	G4LogicalVolume *logicBareSource_Dt, *logicBS_Disk_Dt, *logicBSD_Ring_Dt;
	G4VPhysicalVolume *physBareSource_Dt, *physBS_Disk_Dt, *physBSD_Ring_Dt;
	// End Radioactive Source

	//Liquid Scintillator
	G4Material *matContainer, *matLiquid;
	G4bool isLiquid, isCupDetector;
	G4double container_radius, container_height_half, container_thickness, d_pos_z;
	G4LogicalVolume *logiContainer_F, *logiContainer_B, *logicLiquid_F, *logicLiquid_B;		// F = Front, B = Back
	G4LogicalVolume *logicPlaneDetector_W, *logicPlaneDetector_C, *logicPlaneDetector_L;	// W = World, C = Container, L = Liquid
	G4LogicalVolume *logicRingDetector, *logicTubeDetector;
	G4VPhysicalVolume *physContainer_F, *physContainer_B, *physLiquid_F, *physLiquid_B, *physPlaneDetector, *physRingDetector, *physTubeDetector;
	//End Liquid Scintillator

	//Materials
	G4Material *matScintillator, *matWrapping, *matSiPM;
	G4Material *matXe;  // Xenon gas for test
	// Radioactive Source (Positron Source)
	G4Material *matTi, *matNaCl, *matCsI;
	G4Material *matSi, *matTeflon;



};

#endif
