#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH


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
	// Ideal Detector
	void ConstructShell_Detector();
	// End Ideal Detector
	// tpc
	void ConstructTPC();
	// End tpc
	// Calorimeter
	void ConstructCalorimeter();
	// End Calorimeter
	//Construct source
	void ConstructSource();
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


	//Materials
	G4Material *matXe;  // Xenon gas for test
	// Radioactive Source (Positron Source)
	G4Material *matTi, *matNaCl, *matCsI;


};

#endif
