#include "construction.hh"
#include "CADMesh.hh"
MyDetectorConstruction::MyDetectorConstruction() 
{
	// Define required materials

	DefineMaterials();
	testMaterialName = "LSO";
	isDetector_Shell = true;
	isSource=false;
	isTPC = false;
	isCalorimeter = false;
	// Set the material for each logical volume
	matWorld = Vacuum; //Vacuum;
	// Set the default of each logical volume to be NULL so the sensitive detector selector can work well
	logicDetector_Shell = NULL;
	logicTPC = NULL;
	logicCalorimeter = NULL;

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
	//LSO
	matLSO = new G4Material("Lu2SiO5", 7.4*g/cm3, 3);
	G4Element* Lu = nist->FindOrBuildElement("Lu");
	G4Element* Si = nist->FindOrBuildElement("Si");
	G4Element* O = nist->FindOrBuildElement("O");
	matLSO->AddElement(Lu, 2);
	matLSO->AddElement(Si, 1);
	matLSO->AddElement(O, 5);
	// End LSO
	//LYSO
	matLYSO = new G4Material("Lu2(1-x)Y2xSiO5", 7.2*g/cm3, 4);
	matLYSO->AddElement(Lu, 1);
	G4Element* Y = nist->FindOrBuildElement("Y");
	matLYSO->AddElement(Y, 1);
	matLYSO->AddElement(Si, 1);
	matLYSO->AddElement(O, 5);
	// End LYSO
	//LaBr3
	matLaBr3 = new G4Material("LaBr3", 5.3*g/cm3, 2);
	G4Element* La = nist->FindOrBuildElement("La");
	G4Element* Br = nist->FindOrBuildElement("Br");
	matLaBr3->AddElement(La, 1);
	matLaBr3->AddElement(Br, 3);
	// End LaBr3
	//GAGG
	matGAGG = new G4Material("GAGG", 6.63*g/cm3, 4);
	G4Element* Gd = nist->FindOrBuildElement("Gd");
	G4Element* Al = nist->FindOrBuildElement("Al");
	G4Element* Ga = nist->FindOrBuildElement("Ga");
	matGAGG->AddElement(Gd, 3);
	matGAGG->AddElement(Al, 2);
	matGAGG->AddElement(Ga, 3);
	matGAGG->AddElement(O, 12);
	// End GAGG

	//BGO
	matBGO = new G4Material("BGO", 7.13*g/cm3, 3);
	G4Element* Bi = nist->FindOrBuildElement("Bi");
	G4Element* Ge = nist->FindOrBuildElement("Ge");
	matBGO->AddElement(Bi, 4);
	matBGO->AddElement(Ge, 1);
	matBGO->AddElement(O, 12);
	// End BGO
	//NaI
	matNaI = new G4Material("NaI", 3.67*g/cm3, 2);
	G4Element* Na = nist->FindOrBuildElement("Na");
	G4Element* I = nist->FindOrBuildElement("I");
	matNaI->AddElement(Na, 1);
	matNaI->AddElement(I, 1);	
	// End NaI
	// Ti_
	G4Element* Ti = nist->FindOrBuildElement("Ti");
	matTi = new G4Material("Ti_", 4.507*g/cm3, 1);	//The density of G4_Ti is 4.54*g/cm3
	matTi->AddElement(Ti, 1.);
	//end Ti



	materialMap["LSO"] = matLSO;
    materialMap["LYSO"] = matLYSO;
    materialMap["LaBr3"] = matLaBr3;
    materialMap["GAGG"] = matGAGG;
    materialMap["BGO"] = matBGO;
    materialMap["NaI"] = matNaI;
    materialMap["NaCl"] = matNaCl;
    materialMap["CsI"] = matCsI;
    materialMap["Ti"] = matTi;
    materialMap["Xe"] = matXe;
    materialMap["Air"] = Air;
    materialMap["Vacuum"] = Vacuum;
}
G4Material* MyDetectorConstruction::GetTestMaterial(const G4String& name) {
    auto it = materialMap.find(name);
    if (it != materialMap.end()) {
        return it->second;
    } else {
        G4cerr << "Unknown material name: " << name << ". Defaulting to LSO." << G4endl;
        return matLSO;
    }
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
	fMessenger->DeclareProperty("TestMaterial", testMaterialName, "Set the material name for Shell Detector (valid: LSO, LYSO, LaBr3, GAGG, BGO, NaI, NaCl, CsI, Ti, Xe)");}
// Construct All physical volumes
G4VPhysicalVolume* MyDetectorConstruction::Construct() {
	G4double xWorld = 600*m;
	G4double yWorld = 600*m;
	G4double zWorld = 600*m;

	// A cubic world with volume 1.5 m*1.5 m*1.5 m
	G4Box* solidWorld = new G4Box("solidWorld", xWorld, yWorld, zWorld);
	logicWorld = new G4LogicalVolume(solidWorld, matWorld, "logicWorld");
	physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);
	if (isDetector_Shell)
		ConstructShell_Detector();
	if (isSource)
		ConstructSource();
	if (isTPC)
		ConstructTPC();
	if (isCalorimeter)
		ConstructCalorimeter();

	return physWorld;
}

// Set Sensitive Detector(SD) and Field
void MyDetectorConstruction::ConstructSDandField() {
	G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    Tracker* tracker0 = new Tracker("ShellTracker");
	//Detect_reference *detect_reference = new Detect_reference("Detect_reference");
	sdManager->AddNewDetector(tracker0);
	//sdManager->AddNewDetector(detect_reference);
	if(logicDetector_Shell != NULL)
		logicDetector_Shell->SetSensitiveDetector(tracker0);
	//if(logicBareSource != NULL)
	//	logicBareSource->SetSensitiveDetector(detect_reference);
	//if(logicDisk != NULL)
	//	logicDisk->SetSensitiveDetector(detect_reference);
	//if(logicRing != NULL)
	//	logicRing->SetSensitiveDetector(detect_reference);

}
// Ideal Detector
void MyDetectorConstruction::ConstructShell_Detector() {
	matTest = GetTestMaterial(testMaterialName);
	G4double shell_thickness = 500.*m;
	G4double inner_radius = 0.*cm; //25.*cm+80.*cm;
	G4double outer_radius = inner_radius + shell_thickness;
	G4Sphere* solidDetector_Shell = new G4Sphere("solidDetector_Shell", inner_radius, outer_radius, 0.*deg, 360.*deg, 0.*deg, 360.*deg);
	logicDetector_Shell = new G4LogicalVolume(solidDetector_Shell, matTest, "logicDetector_Shell");
	physDetector_Shell = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, 0.*m), logicDetector_Shell, "Detector_Shell", logicWorld, false, 0, true);
}
// End Ideal Detector
void MyDetectorConstruction::ConstructTPC() {
	G4double inner_radius = 0.*m;
	G4double outer_radius = 25.*cm;
	G4Sphere* solidDetector_Shell = new G4Sphere("solidTPC", inner_radius, outer_radius, 0.*deg, 360.*deg, 0.*deg, 360.*deg);
	logicTPC = new G4LogicalVolume(solidDetector_Shell, matXe, "logicTPC");
	physTPC = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, 0.*m), logicTPC, "DetectorTPC", logicWorld, false, 0, true);
}
// End Ideal Detector
void MyDetectorConstruction::ConstructCalorimeter() {
	G4double inner_radius = 25.*cm;
	G4double outer_radius = 25.*cm+80.*cm;
	G4Sphere* solidDetector_Shell = new G4Sphere("solidCalorimeter", inner_radius, outer_radius, 0.*deg, 360.*deg, 0.*deg, 360.*deg);
	logicCalorimeter = new G4LogicalVolume(solidDetector_Shell, matCsI, "logicCalorimeter");
	physCalorimeter = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, 0.*m), logicCalorimeter, "Calorimeter", logicWorld, false, 0, true);
}
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