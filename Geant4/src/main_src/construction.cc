#include "construction.hh"
#include "CADMesh.hh"
MyDetectorConstruction::MyDetectorConstruction() {
	// Define required materials

	DefineMaterials();

	isDetector_Shell = true;
	isSource=false;
	// Set the material for each logical volume
	matWorld = Vacuum; //Vacuum;

	// Set the default of each logical volume to be NULL so the sensitive detector selector can work well
	logicDetector_Shell = NULL;

	DefineMessenger();
}
MyDetectorConstruction::~MyDetectorConstruction()
{}

G4String MyDetectorConstruction::file_name = "";

void MyDetectorConstruction::DefineMaterials() {
	G4NistManager* nist = G4NistManager::Instance();
	
	Air = nist->FindOrBuildMaterial("G4_AIR");
	Vacuum = nist->FindOrBuildMaterial("G4_Galactic");
	// Defining Xenon gas for test
	auto a = 131.29*g/mole;
	G4Element* Xe = new G4Element("Xe", "Xe", 54., a);
	auto density = 5.858*mg/cm3;  
	double pressure = 1*bar;  // [X->Your choice]
	double temperature = 296.15*kelvin;  // [your choice]
	matXe  = new G4Material("matXe", density, 1, kStateGas, temperature, pressure);
	matXe->AddElement(Xe, 1);  //--> Monoatomic nature

	//NaCl
	matNaCl = new G4Material("NaCl", 2.16*g/cm3, 2);
	matNaCl->AddElement(nist->FindOrBuildElement("Na"), 1);
	matNaCl->AddElement(nist->FindOrBuildElement("Cl"), 4);
	// End NaCl
	// CsI
	matCsI = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");
	// End CsI

	// Ti_
	G4Element* Ti = nist->FindOrBuildElement("Ti");
	matTi = new G4Material("Ti_", 4.507*g/cm3, 1);	//The density of G4_Ti is 4.54*g/cm3
	matTi->AddElement(Ti, 1.);
	//end Ti
}

void MyDetectorConstruction::DefineMessenger() {
	// The A is a placeholder for the user defined commands (fMessenger)
	G4int placeHolder = 0;
	// These are user defined commands for use in User-Interface(UI) mode and batch mode(using macro file)
	fMessenger = new G4GenericMessenger(this, "/MyDetector/", "Macros");
	fMessenger->DeclareProperty("control/execute region_setup.mac", placeHolder, "Set the active region (cylinder locate at origin, radius = 9.53/2*mm, half height = 0.0001*mm)");
	fMessenger->DeclareProperty("control/execute rebuild.mac",placeHolder,"Rebuild Selected Physical Volume inside a 1.5*1.5*1.5 m^3 Cubic World contains Air, its center is the origin");
	fMessenger->DeclareProperty("isDetector_Shell", isDetector_Shell, "Construct Shell Detector (spherical shell locate at origin, inner radius = 3*cm, thickness = 1*nm)");
	fMessenger->DeclareProperty("setFileName", file_name, "Set the name of output root file");
}
// Construct All physical volumes
G4VPhysicalVolume* MyDetectorConstruction::Construct() {
	G4double xWorld = 2*m;
	G4double yWorld = 2*m;
	G4double zWorld = 2*m;

	// A cubic world with volume 1.5 m*1.5 m*1.5 m
	G4Box* solidWorld = new G4Box("solidWorld", xWorld, yWorld, zWorld);
	logicWorld = new G4LogicalVolume(solidWorld, matWorld, "logicWorld");
	physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);
	if (isDetector_Shell)
		ConstructShell_Detector();
	if (isSource)
		ConstructSource();
	return physWorld;
}

// Set Sensitive Detector(SD) and Field
void MyDetectorConstruction::ConstructSDandField() {
	G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    Tracker* tracker0 = new Tracker("ShellTracker");
	Detect_reference *detect_reference = new Detect_reference("Detect_reference");
	sdManager->AddNewDetector(tracker0);
	sdManager->AddNewDetector(detect_reference);
	if(logicDetector_Shell != NULL)
		logicDetector_Shell->SetSensitiveDetector(tracker0);
	if(logicBareSource != NULL)
		logicBareSource->SetSensitiveDetector(detect_reference);
	if(logicDisk != NULL)
		logicDisk->SetSensitiveDetector(detect_reference);
	if(logicRing != NULL)
		logicRing->SetSensitiveDetector(detect_reference);

}
// Ideal Detector
void MyDetectorConstruction::ConstructShell_Detector() {
	G4double shell_thickness = 1.*m;
	G4double inner_radius = 0.*cm;
	G4double outer_radius = inner_radius + shell_thickness;
	G4Sphere* solidDetector_Shell = new G4Sphere("solidDetector_Shell", inner_radius, outer_radius, 0.*deg, 360.*deg, 0.*deg, 360.*deg);
	logicDetector_Shell = new G4LogicalVolume(solidDetector_Shell, matXe, "logicDetector_Shell");
	physDetector_Shell = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, 0.*m), logicDetector_Shell, "Detector_Shell", logicWorld, false, 0, true);
}
// End Ideal Detector

//Construct source
void MyDetectorConstruction::ConstructSource(){
	ring_radius = 19.1/2*mm;
	ring_height_half = 0.254*mm;			//Supported by two 0.254 mm Ti disks
	disk_radius = 9.53/2*mm;
	disk_height_half = 0.00508*mm;			//The activity is placed between two layers of 0.00508*mm Ti foil which is 0.0102*mm in total
	bare_source_radius = disk_radius;
	bare_source_height_half = 0.0001*mm;	//The thickness of bare source is not provide so this is a made up value
	G4Tubs* solidRing = new G4Tubs("solidRing", disk_radius, ring_radius, ring_height_half, 0.*deg, 360.*deg);
	logicRing = new G4LogicalVolume(solidRing, matTi, "logicRing");
	G4Translate3D trans(G4ThreeVector(0*m, 0*m, 0*cm));
	G4Rotate3D rotY(90.*deg, G4ThreeVector(0.,1.,0.));
	G4Transform3D tranTest = rotY*trans;
	physRing = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, 0.*m), logicRing, "Ring", logicWorld, false, 0, true);

	G4Tubs* solidDisk = new G4Tubs("solidDisk", 0.*nm, disk_radius, disk_height_half, 0.*deg, 360.*deg);
	logicDisk = new G4LogicalVolume(solidDisk, matTi, "logicDisk");
	physDisk = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, 0.*m), logicDisk, "Disk", logicWorld, false, 0, true);

	G4Tubs* solidBareSource = new G4Tubs("solidBareSource", 0.*nm, bare_source_radius, bare_source_height_half, 0.*deg, 360.*deg);
	logicBareSource = new G4LogicalVolume(solidBareSource, matNaCl, "logicBareSource");
	physBareSource = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, 0.*m), logicBareSource, "BareSource", logicWorld, false, 0, true);

}
// End Construct source