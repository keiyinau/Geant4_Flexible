#include "construction.hh"
#include "CADMesh.hh"
MyDetectorConstruction::MyDetectorConstruction() {
	// Define required materials

	DefineMaterials();

	isDetector_Shell = false;
	isSource=false;
	isTPC = false;
	isCalorimeter = true;
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
bool MyDetectorConstruction::readAndProcessData(const std::string& filename, 
                       std::vector<double>& emission_Energy, 
                       std::vector<double>& emission_fractions) {
    // Open the file
    std::ifstream datafile(filename);
    if (!datafile.is_open()) {
        return false; // Return false if file cannot be opened
    }

    // Read and parse the file
    std::string line;
    while (std::getline(datafile, line)) {
        std::istringstream iss(line);
        double wlen, fraction;
        char delim;
        if (iss >> wlen >> delim >> fraction && delim == ',') {
            emission_Energy.push_back(1239.84193 / wlen); // E=hc/λ
            emission_fractions.push_back(fraction);
        }
    }

    // If no data was read, return false
    if (emission_Energy.empty()) {
        return false;
    }

    // Pair energies and fractions for sorting
    std::vector<std::pair<double, double>> paired;
    for (size_t i = 0; i < emission_Energy.size(); ++i) {
        paired.emplace_back(emission_Energy[i], emission_fractions[i]);
    }

    // Sort by energy (first element)
    std::sort(paired.begin(), paired.end());

    // Update vectors with sorted values
    emission_Energy.clear();
    emission_fractions.clear();
    for (const auto& p : paired) {
        emission_Energy.push_back(p.first);
        emission_fractions.push_back(p.second);
    }
	

    return true; // Success
	}
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
	G4MaterialPropertiesTable* mptCsI = new G4MaterialPropertiesTable();
	//CsI Emssion spectrum
	std::vector<double> CsI_emission_Energy, CsI_emission_fractions;
	readAndProcessData("EmissionSpectrum_295K.csv", CsI_emission_Energy, CsI_emission_fractions);
	mptCsI->AddConstProperty("RESOLUTIONSCALE", 1.);	
	mptCsI->AddProperty("SCINTILLATIONCOMPONENT1", CsI_emission_Energy, CsI_emission_fractions);	
	mptCsI->AddConstProperty("SCINTILLATIONYIELD", 30./keV);
	mptCsI->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 20.0*ns);	
	matCsI->SetMaterialPropertiesTable(mptCsI);
	// End CsI

	// Define Tapflon(teflon) for wrapping
	matTeflon = nist->FindOrBuildMaterial("G4_TEFLON");
	std::vector<double> tapflon_reflectance_Energy, tapflon_reflectance_fractions;
	readAndProcessData("EmissionSpectrum_295K.csv", tapflon_reflectance_Energy, tapflon_reflectance_fractions);
	G4MaterialPropertiesTable* matTeflon = new G4MaterialPropertiesTable();
	matTeflon->AddProperty("REFLECTIVITY", tapflon_reflectance_Energy, tapflon_reflectance_fractions);
	// End Tapflon


	// Define SiPM
	matSi = nist->FindOrBuildMaterial("G4_Si");

	// End SiPM



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
	Calorimeter* calorimeter = new Calorimeter("Calorimeter");
	//Detect_reference *detect_reference = new Detect_reference("Detect_reference");
	sdManager->AddNewDetector(tracker0);
	sdManager->AddNewDetector(calorimeter);
	//sdManager->AddNewDetector(detect_reference);
	if(logicDetector_Shell != NULL)
		logicDetector_Shell->SetSensitiveDetector(tracker0);
	if(logicCalorimeter!=NULL)
		logicCalorimeter->SetSensitiveDetector(calorimeter);
	//if(logicBareSource != NULL)
	//	logicBareSource->SetSensitiveDetector(detect_reference);
	//if(logicDisk != NULL)
	//	logicDisk->SetSensitiveDetector(detect_reference);
	//if(logicRing != NULL)
	//	logicRing->SetSensitiveDetector(detect_reference);

}
// Ideal Detector
void MyDetectorConstruction::ConstructShell_Detector() {
	G4double shell_thickness = 1.*nm;
	G4double inner_radius = 25.*cm+80.*cm;
	G4double outer_radius = inner_radius + shell_thickness;
	G4Sphere* solidDetector_Shell = new G4Sphere("solidDetector_Shell", inner_radius, outer_radius, 0.*deg, 360.*deg, 0.*deg, 360.*deg);
	logicDetector_Shell = new G4LogicalVolume(solidDetector_Shell, matCsI, "logicDetector_Shell");
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
	std::string Scintillator_name_list[] = {"CsI_CenteredCsI_CenteredS1(1)"};
    std::string SiPM_name_list[] = {"CsI_CenteredCsI_CenteredD1(1)"};
    std::string Tapflon_name_list[] = {"CsI_CenteredCsI_CenteredT1(1)"};
    int Size_of_Scintillator_name_list = sizeof(Scintillator_name_list)/sizeof(std::string);
    int Size_of_SiPM_name_list = sizeof(SiPM_name_list)/sizeof(std::string);
    int Size_of_Tapflon_name_list = sizeof(Tapflon_name_list)/sizeof(std::string);
    std::vector<G4LogicalVolume*> logicScintillators(Size_of_Scintillator_name_list);
	std::vector<G4LogicalVolume*> logicTapflon(Size_of_Tapflon_name_list);
    std::vector<G4LogicalVolume*> logicSiPM(Size_of_SiPM_name_list);
    std::vector<G4VPhysicalVolume*> physScintillators(Size_of_Scintillator_name_list);
    std::vector<G4VPhysicalVolume*> physTapflon(Size_of_Tapflon_name_list);
    std::vector<G4VPhysicalVolume*> physSiPM(Size_of_SiPM_name_list);
    G4double angle = 90 * deg;
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angle);
	for (int i = 0; i < Size_of_Scintillator_name_list; i++) {
        std::string name_scint = Scintillator_name_list[i];
        auto scintillatorDet = CADMesh::TessellatedMesh::FromSTL(name_scint + ".stl");
        auto ScintillatorDet = scintillatorDet->GetSolid();
        G4LogicalVolume* logicScintillator_pre = new G4LogicalVolume(ScintillatorDet, matCsI, name_scint + "Logic");
        logicScintillators[i] = logicScintillator_pre;
        physScintillators[i] = new G4PVPlacement(rotation, G4ThreeVector(), logicScintillator_pre, name_scint, logicWorld, false, i, true);    
    }
	for (int i = 0; i < Size_of_SiPM_name_list; i++) {
        std::string name_scint = SiPM_name_list[i];
        auto scintillatorDet = CADMesh::TessellatedMesh::FromSTL(name_scint + ".stl");
        auto ScintillatorDet = scintillatorDet->GetSolid();
        G4LogicalVolume* logicSiPM_pre = new G4LogicalVolume(ScintillatorDet, matCsI, name_scint + "Logic");
		logicCalorimeter=logicSiPM_pre;
        logicSiPM[i] = logicCalorimeter;
        physSiPM[i] = new G4PVPlacement(rotation, G4ThreeVector(), logicCalorimeter, name_scint, logicWorld, false, i, true);    
    }
	for (int i = 0; i < Size_of_Tapflon_name_list; i++) {
        std::string name_scint = Tapflon_name_list[i];
        auto scintillatorDet = CADMesh::TessellatedMesh::FromSTL(name_scint + ".stl");
        auto ScintillatorDet = scintillatorDet->GetSolid();
        G4LogicalVolume* logicTapflon_pre = new G4LogicalVolume(ScintillatorDet, matTeflon, name_scint + "Logic");
        logicTapflon[i] = logicTapflon_pre;
        physTapflon[i] = new G4PVPlacement(rotation, G4ThreeVector(), logicTapflon_pre, name_scint, logicWorld, false, i, true);    
    }

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