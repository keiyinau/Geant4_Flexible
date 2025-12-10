#include "construction.hh"
#include "CADMesh.hh"
MyDetectorConstruction::MyDetectorConstruction() {
	// Define required materials

	DefineMaterials();

	isDetector_Shell = false;
	isSource=true;
	isTPC = false;
	isCalorimeter = true;
    isLiquid=true;
    is3DCalorimeter=true;
	// Set the material for each logical volume
	matWorld = Air; //Vacuum;
    matLiquid=matWater;
    matContainer=matAcrylic;
    matScintillator=matLSO;
    matSiPM=matSi;
    matWrapping=matTeflon;
	// Set the default of each logical volume to be NULL so the sensitive detector selector can work well
	logicDetector_Shell = NULL;
	logicTPC = NULL;
	logicCalorimeter = NULL;


    ring_radius = 19.1/2*mm;
	ring_height_half = 0.254*mm;			//Supported by two 0.254 mm Ti disks
	disk_radius = 9.53/2*mm;
	disk_height_half = 0.00508*mm;			//The activity is placed between two layers of 0.00508*mm Ti foil which is 0.0102*mm in total
	bare_source_radius = disk_radius;
	bare_source_height_half = 0.0001*mm;	//The thickness of bare source is not provide so this is a made up value
    container_radius = 2*cm;				//All variables of container are made up.
	container_height_half = 2/2*cm;
	container_thickness = 1.*mm;
	d_pos_z = 0.02*mm;						//Distance between two nearest plane detectors, spacing of two plane detectors

	DefineMessenger();




}
MyDetectorConstruction::~MyDetectorConstruction()
{}
bool MyDetectorConstruction::readAndProcessData(const std::string& filename, 
                       std::vector<G4double>& emission_Energy, 
                       std::vector<G4double>& emission_fractions) {
    // Open the file
    std::ifstream datafile(filename);
    if (!datafile.is_open()) {
        return false; // Return false if file cannot be opened
    }

    // Read and parse the file
    std::string line;
    while (std::getline(datafile, line)) {
        std::istringstream iss(line);
        G4double wlen, fraction;
        char delim;
        if (iss >> wlen >> delim >> fraction && delim == ',') {
            emission_Energy.push_back(wlen); // E=hc/λ
            emission_fractions.push_back(fraction);
        }
    }

    // If no data was read, return false
    if (emission_Energy.empty()) {
        return false;
    }

    // Pair energies and fractions for sorting
    std::vector<std::pair<G4double, G4double>> paired;
    for (size_t i = 0; i < emission_Energy.size(); ++i) {
        paired.emplace_back(emission_Energy[i], emission_fractions[i]);
    }

    // Sort by energy (first element)
    std::sort(paired.begin(), paired.end());

    // Update vectors with sorted values
    emission_Energy.clear();
    emission_fractions.clear();
    for (const auto& p : paired) {
        emission_Energy.push_back(p.first*eV);
        emission_fractions.push_back(p.second);
    }
	

    return true; // Success
	}

bool MyDetectorConstruction::readAndProcessData_Energy(const std::string& filename, 
                       std::vector<G4double>& emission_Energy, 
                       std::vector<G4double>& emission_fractions) {
    // Open the file
    std::ifstream datafile(filename);
    if (!datafile.is_open()) {
        return false; // Return false if file cannot be opened
    }

    // Read and parse the file
    std::string line;
    while (std::getline(datafile, line)) {
        std::istringstream iss(line);
        G4double wlen, fraction;
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
    std::vector<std::pair<G4double, G4double>> paired;
    for (size_t i = 0; i < emission_Energy.size(); ++i) {
        paired.emplace_back(emission_Energy[i], emission_fractions[i]);
    }

    // Sort by energy (first element)
    std::sort(paired.begin(), paired.end());

    // Update vectors with sorted values
    emission_Energy.clear();
    emission_fractions.clear();
    for (const auto& p : paired) {
        emission_Energy.push_back(p.first*eV);
        emission_fractions.push_back(p.second);
    }
	

    return true; // Success
	}

bool MyDetectorConstruction::readAndProcessData_txt(const std::string& filename, 
                                               std::vector<G4double>& emission_Energy, 
                                               std::vector<G4double>& emission_fractions) {
   std::ifstream datafile(filename);
    if (!datafile) {
        std::cerr << "Error: Cannot open file " << filename << "\n";
        return false;
    }

    std::vector<std::pair<G4double, G4double>> paired;
    std::string line;
    while (std::getline(datafile, line)) {
        std::istringstream iss(line);
        G4double wlen, fraction;
        if (iss >> wlen >> fraction) { // Space-separated values
            if (wlen <= 0) continue; // Skip invalid wavelengths
            paired.emplace_back(wlen, fraction); // Energy (eV), fraction
        }
    }

    if (paired.empty()) {
        std::cerr << "Error: No valid data read from " << filename << "\n";
        return false;
    }

    std::sort(paired.begin(), paired.end()); // Sort by energy (increasing)

    emission_Energy.clear();
    emission_fractions.clear();
    for (const auto& p : paired) {
        emission_Energy.push_back(p.first*eV);
        emission_fractions.push_back(p.second);
    }

    return true;
}
bool MyDetectorConstruction::readAndProcessData_Energy_txt(const std::string& filename, 
                                                      std::vector<G4double>& emission_Energy, 
                                                      std::vector<G4double>& emission_fractions) {
    std::ifstream datafile(filename);
    if (!datafile) {
        std::cerr << "Error: Cannot open file " << filename << "\n";
        return false;
    }

    std::vector<std::pair<G4double, G4double>> paired;
    std::string line;
    while (std::getline(datafile, line)) {
        std::istringstream iss(line);
        G4double wlen, fraction;
        if (iss >> wlen >> fraction) { // Space-separated values
            if (wlen <= 0) continue; // Skip invalid wavelengths
            paired.emplace_back(1239.84193 / wlen, fraction); // Energy (eV), fraction
        }
    }

    if (paired.empty()) {
        std::cerr << "Error: No valid data read from " << filename << "\n";
        return false;
    }

    std::sort(paired.begin(), paired.end()); // Sort by energy (increasing)

    emission_Energy.clear();
    emission_fractions.clear();
    for (const auto& p : paired) {
        emission_Energy.push_back(p.first*eV);
        emission_fractions.push_back(p.second);
    }

    return true;
}
bool MyDetectorConstruction::readAndProcessData_Energy_cm_txt(const std::string& filename, 
                                                      std::vector<G4double>& emission_Energy, 
                                                      std::vector<G4double>& emission_fractions) {
    std::ifstream datafile(filename);
    if (!datafile) {
        std::cerr << "Error: Cannot open file " << filename << "\n";
        return false;
    }

    std::vector<std::pair<G4double, G4double>> paired;
    std::string line;
    while (std::getline(datafile, line)) {
        std::istringstream iss(line);
        G4double wlen, fraction;
        if (iss >> wlen >> fraction) { // Space-separated values
            if (wlen <= 0) continue; // Skip invalid wavelengths
            paired.emplace_back(1239.84193 / wlen, fraction); // Energy (eV), fraction
        }
    }

    if (paired.empty()) {
        std::cerr << "Error: No valid data read from " << filename << "\n";
        return false;
    }

    std::sort(paired.begin(), paired.end()); // Sort by energy (increasing)

    emission_Energy.clear();
    emission_fractions.clear();
    for (const auto& p : paired) {
        emission_Energy.push_back(p.first*eV);
        emission_fractions.push_back(p.second*cm);
    }

    return true;
}

G4String MyDetectorConstruction::file_name = "";

void MyDetectorConstruction::DefineMaterials() {
	G4NistManager* nist = G4NistManager::Instance();
	//Define the world material as Air
	Air = nist->FindOrBuildMaterial("G4_AIR");
    std::vector<G4double> Air_absorption_Energy, Air_absorption_Index;
    readAndProcessData_Energy_cm_txt("AbsorptionLength_Air.txt", Air_absorption_Energy, Air_absorption_Index);
    G4MaterialPropertiesTable* mptAir = new G4MaterialPropertiesTable();
    mptAir->AddProperty("RINDEX", "Air");
    mptAir->AddProperty("ABSLENGTH", Air_absorption_Energy, Air_absorption_Index,Air_absorption_Index.size());
    //Air->SetMaterialPropertiesTable(mptAir);

    // Define the world material as vacuum
	Vacuum = nist->FindOrBuildMaterial("G4_Galactic");
	// Defining Xenon gas for test
	auto a = 131.29*g/mole;
	G4Element* Xe = new G4Element("Xe", "Xe", 54., a);
	auto density = 5.858*mg/cm3;  
	double pressure = 1*bar;  // [X->Your choice]
	double temperature = 296.15*kelvin;  // [your choice]
	matXe  = new G4Material("matXe", density, 1, kStateGas, temperature, pressure);
	matXe->AddElement(Xe, 1);  //--> Monoatomic nature
    // End Xenon gas

    // Define water
    matWater = nist->FindOrBuildMaterial("G4_WATER");
    G4MaterialPropertiesTable* mptWater = new G4MaterialPropertiesTable();
    std::vector<G4double> Water_absorption_Energy, Water_absorption_Index;
    readAndProcessData_Energy_cm_txt("AbsorptionLength_Water.txt", Water_absorption_Energy, Water_absorption_Index);
    mptWater->AddProperty("RINDEX", "Water");
    mptWater->AddProperty("ABSLENGTH", Water_absorption_Energy, Water_absorption_Index,Water_absorption_Index.size());
    //matWater->SetMaterialPropertiesTable(mptWater);
    // End water



	//NaCl
	matNaCl = new G4Material("NaCl", 2.16*g/cm3, 2);
	matNaCl->AddElement(nist->FindOrBuildElement("Na"), 1);
	matNaCl->AddElement(nist->FindOrBuildElement("Cl"), 4);
	// End NaCl
    matLSO = new G4Material("Lu2SiO5", 7.4*g/cm3, 3);
	G4Element* Lu = nist->FindOrBuildElement("Lu");
	G4Element* Si = nist->FindOrBuildElement("Si");
	G4Element* O = nist->FindOrBuildElement("O");
	matLSO->AddElement(Lu, 2);
	matLSO->AddElement(Si, 1);
	matLSO->AddElement(O, 5);
    G4MaterialPropertiesTable* mptLSO = new G4MaterialPropertiesTable();
    //LSO Emssion spectrum
    std::vector<G4double> LSO_emission_Energy, LSO_emission_fractions;
    readAndProcessData_Energy("EmissionSpectrum_LSO_Ce_295K.csv", LSO_emission_Energy, LSO_emission_fractions);
    std::vector<G4double> LSO_refraction_Energy, LSO_refraction_Index;
    readAndProcessData_txt("RefractiveIndex_LSO_Ce.txt", LSO_refraction_Energy, LSO_refraction_Index);
    std::vector<G4double> LSO_absorption_Energy, LSO_absorption_Index;
    readAndProcessData_Energy_cm_txt("AbsorptionLength_LSO_Ce.txt", LSO_absorption_Energy, LSO_absorption_Index);
    mptLSO->AddConstProperty("RESOLUTIONSCALE", 0.);
    mptLSO->AddProperty("SCINTILLATIONCOMPONENT1", LSO_emission_Energy, LSO_emission_fractions,LSO_emission_fractions.size());
    mptLSO->AddProperty("RINDEX", LSO_refraction_Energy, LSO_refraction_Index,LSO_refraction_Index.size());
    mptLSO->AddProperty("ABSLENGTH", LSO_absorption_Energy, LSO_absorption_Index,LSO_absorption_Index.size());
    mptLSO->AddConstProperty("SCINTILLATIONYIELD", 26/keV);
    mptLSO->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 40.0*ns);
    //matLSO->SetMaterialPropertiesTable(mptLSO);
    // End LSO



	// CsI
	matCsI = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");
	G4MaterialPropertiesTable* mptCsI = new G4MaterialPropertiesTable();
	//CsI Emssion spectrum
	std::vector<G4double> CsI_emission_Energy, CsI_emission_fractions;
	readAndProcessData_Energy("EmissionSpectrum_295K.csv", CsI_emission_Energy, CsI_emission_fractions);
	std::vector<G4double> CsI_refraction_Energy, CsI_refraction_Index;
	readAndProcessData_txt("RefractiveIndexINFO_CsI.txt", CsI_refraction_Energy, CsI_refraction_Index);
	std::vector<G4double> CsI_transmission_Energy, CsI_rtransmission_Index;
	readAndProcessData_txt("transmittance_CsI.txt", CsI_transmission_Energy, CsI_rtransmission_Index);
	std::vector<G4double> CsI_absorption_Energy, CsI_absorption_Index;
	readAndProcessData_Energy_cm_txt("Absorption_CsITi.txt", CsI_absorption_Energy, CsI_absorption_Index);



	mptCsI->AddConstProperty("RESOLUTIONSCALE", 1.);	
	mptCsI->AddProperty("SCINTILLATIONCOMPONENT1", CsI_emission_Energy, CsI_emission_fractions,CsI_emission_fractions.size());
	mptCsI->AddProperty("RINDEX", CsI_refraction_Energy, CsI_refraction_Index,CsI_refraction_Index.size());	
	mptCsI->AddProperty("TRANSMITTANCE", CsI_transmission_Energy, CsI_rtransmission_Index,CsI_rtransmission_Index.size());	
    mptCsI->AddProperty("ABSLENGTH", CsI_absorption_Energy, CsI_absorption_Index,CsI_absorption_Index.size());
	mptCsI->AddConstProperty("SCINTILLATIONYIELD", 3./keV);
	mptCsI->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 25.0*ns);	
	//matCsI->SetMaterialPropertiesTable(mptCsI);
	// End CsI

	// Define Aluminium for wrapping and protection
	matAl = nist->FindOrBuildMaterial("G4_Al");
	G4MaterialPropertiesTable* mptAl = new G4MaterialPropertiesTable();
    const G4int nEntries = 2; // Example with two points
    G4double PhotonEnergy[nEntries] = {1.5 * eV, 3.0 * eV}; // Example energy range
    G4double RIndex_al[nEntries] = {1.37, 0.44}; // Example refractive index values
    // Add the properties to the table
    mptAl->AddProperty("RINDEX", PhotonEnergy, RIndex_al, nEntries);
    //matAl->SetMaterialPropertiesTable(mptAl);
    // End Aluminium

	// Define Acrylic
	matAcrylic = nist->FindOrBuildMaterial("G4_PLEXIGLASS");
	G4MaterialPropertiesTable* mptAcrylic = new G4MaterialPropertiesTable();
    mptAcrylic->AddProperty("RINDEX", "PMMA");
	//matAcrylic->SetMaterialPropertiesTable(mptAcrylic);
	// End Acrylic



	// Define Tapflon(teflon) for wrapping
	matTeflon = nist->FindOrBuildMaterial("G4_TEFLON");
	std::vector<G4double> tapflon_reflectance_Energy, tapflon_reflectance_fractions;
	readAndProcessData_Energy_txt("teflon_Reflectance-modified.txt", tapflon_reflectance_Energy, tapflon_reflectance_fractions);
    std::vector<G4double> tapflon_refraction_Energy, tapflon_refraction_Index;
	readAndProcessData_Energy_txt("Refraction_Index_Teflon_Gray.txt", tapflon_refraction_Energy, tapflon_refraction_Index);
	G4MaterialPropertiesTable* mptTeflon = new G4MaterialPropertiesTable();
	mptTeflon->AddProperty("REFLECTIVITY", tapflon_reflectance_Energy, tapflon_reflectance_fractions,tapflon_reflectance_fractions.size());
    mptTeflon->AddProperty("RINDEX", tapflon_refraction_Energy, tapflon_refraction_Index,tapflon_refraction_Index.size());
    //matTeflon->SetMaterialPropertiesTable(mptTeflon);
    // End Tapflon

	// Define SiPM
	matSi = nist->FindOrBuildMaterial("G4_Si");
	std::vector<G4double> Si_reflectance_Energy, Si_reflectance_fractions;
	readAndProcessData_Energy_txt("Reflectance_Si.txt", Si_reflectance_Energy, Si_reflectance_fractions);
	std::vector<G4double> Si_transmission_Energy, Si_rtransmission_Index;
	readAndProcessData_txt("transmittance_Si.txt", Si_transmission_Energy, Si_rtransmission_Index);
	std::vector<G4double> Si_refraction_Energy, Si_refraction_Index;
	readAndProcessData_txt("RefractiveIndexINFO_Si.txt", Si_refraction_Energy, Si_refraction_Index);

	G4MaterialPropertiesTable* mptSi = new G4MaterialPropertiesTable();
	mptSi->AddProperty("REFLECTIVITY", Si_reflectance_Energy, Si_reflectance_fractions,Si_reflectance_fractions.size());
    mptSi->AddProperty("TRANSMITTANCE", Si_transmission_Energy, Si_rtransmission_Index,Si_rtransmission_Index.size());	
    mptSi->AddProperty("RINDEX", Si_refraction_Energy, Si_refraction_Index,Si_refraction_Energy.size());	
	//matSi->SetMaterialPropertiesTable(mptSi);

    // CsI-Teflon (reflective surface)
    surfCsI_Teflon = new G4OpticalSurface("CsI_Teflon_Surface");
    surfCsI_Teflon->SetType(dielectric_metal); // Teflon as reflective surface
    surfCsI_Teflon->SetModel(unified);
    surfCsI_Teflon->SetFinish(polished);
    //End surface

    // CsI-SiPM (dielectric-dielectric interface)
    surfCsI_SiPM = new G4OpticalSurface("CsI_SiPM_Surface");
    surfCsI_SiPM->SetType(dielectric_dielectric);
    surfCsI_SiPM->SetModel(glisur); // Glisur for smooth dielectric interface
    surfCsI_SiPM->SetFinish(polished);
	// End SiPM

    // CsI-AlFoil (reflective surface)
    surfCsI_AlFoil = new G4OpticalSurface("CsI_AlFoil_Surface");
    surfCsI_AlFoil->SetType(dielectric_metal); // Al as reflective surface
    surfCsI_AlFoil->SetModel(unified);
    surfCsI_AlFoil->SetFinish(ground);
    //End surface



    std::cout<<"==========================="<<std::endl;
    std::cout<<"Printing the material properties of CsI"<<std::endl;
    mptCsI->DumpTable();
    std::cout<<"==========================="<<std::endl;
    std::cout<<"==========================="<<std::endl;
    std::cout<<"Printing the material properties of Al"<<std::endl;
    mptAl->DumpTable();
    std::cout<<"==========================="<<std::endl;
    std::cout<<"Printing the material properties of Si"<<std::endl;
    mptSi->DumpTable();
    std::cout<<"==========================="<<std::endl;



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
	G4double xWorld = 0.22*m;
	G4double yWorld = 0.22*m;
	G4double zWorld = 0.22*m;

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
    if (isLiquid)
        ConstructLiquidScintillator();
    std::cout<<"==========================="<<std::endl;
    std::cout<<"Test if there are overlap"<<std::endl;
    physWorld->CheckOverlaps();
    std::cout<<"==========================="<<std::endl;

	return physWorld;
}

// Set Sensitive Detector(SD) and Field
void MyDetectorConstruction::ConstructSDandField() {
	G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    Tracker* tracker0 = new Tracker("Detector_Real");
    Detect_reference* detect_reference = new Detect_reference("Detect_reference");
	Calorimeter* calorimeter = new Calorimeter("Calorimeter");
    Detect_edep* detect_edep = new Detect_edep("Detector_edep");
	sdManager->AddNewDetector(tracker0);
	sdManager->AddNewDetector(calorimeter);
    sdManager->AddNewDetector(detect_reference);
	sdManager->AddNewDetector(detect_edep);
	if(logicDetector_Shell != NULL)
		logicDetector_Shell->SetSensitiveDetector(detect_reference);
	if(logicCalorimeter!=NULL)
        for(int i=0; i < logicSiPM.size(); i++) {
            logicSiPM[i]->SetSensitiveDetector(calorimeter);
        }
        for(int i=0; i < logicScintillators.size(); i++) {
            logicScintillators[i]->SetSensitiveDetector(detect_edep);
        }
		//logicCalorimeter->SetSensitiveDetector(calorimeter);
	//if(logicBareSource != NULL)
	//	logicBareSource->SetSensitiveDetector(detect_reference);
	//if(logicDisk != NULL)
	//	logicDisk->SetSensitiveDetector(detect_reference);
	//if(logicRing != NULL)
	//	logicRing->SetSensitiveDetector(detect_reference);

}
// Ideal Detector
void MyDetectorConstruction::ConstructShell_Detector() {
	G4double shell_thickness = 1.*nm;//1.*nm;
	G4double inner_radius =50.0*cm;// 25.*cm+80.*cm;
	G4double outer_radius = inner_radius + shell_thickness;
	G4Sphere* solidDetector_Shell = new G4Sphere("solidDetector_Shell", inner_radius, outer_radius, 0.*deg, 360.*deg, 0.*deg, 360.*deg);
	logicDetector_Shell = new G4LogicalVolume(solidDetector_Shell, matWorld, "logicDetector_Shell");
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
void MyDetectorConstruction::ConstructCalorimeter_unit(G4ThreeVector translation,
                                                       G4double angle,
                                                       G4String name)
{
    // --------------------------------------------------------------
    // 1. Rotation
    // --------------------------------------------------------------
    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateX(angle);

    // --------------------------------------------------------------
    // 2. Name lists
    // --------------------------------------------------------------
    const std::string Scintillator_name_list[] = {"CsI"};
    const std::string SiPM_name_list[]       = {"SiPM0"};//,"SiPM1"};//,"SiPM2","SiPM3"};
    const std::string Tapflon_name_list[]    = {"AlFoil"};
    const std::string Protection_name_list[] = {};//{"Bottom","Left","Right","Top"};
    const std::string Acrylic_name_list[]    = {"Acrylic"};

    const int nScint = sizeof(Scintillator_name_list)/sizeof(std::string);
    const int nSiPM  = sizeof(SiPM_name_list)/sizeof(std::string);
    const int nFoil  = sizeof(Tapflon_name_list)/sizeof(std::string);
    const int nProt  = sizeof(Protection_name_list)/sizeof(std::string);
    const int nAcry  = sizeof(Acrylic_name_list)/sizeof(std::string);

    std::vector<G4VPhysicalVolume*> physScint(nScint), physFoil(nFoil),
                                    physSiPM(nSiPM),   physProt(nProt),
                                    physAcry(nAcry);

    // --------------------------------------------------------------
    // 3. CsI crystal (5×5×50 cm³)
    // --------------------------------------------------------------
    for (int i=0;i<nScint;i++) {
        G4Box* box = new G4Box(Scintillator_name_list[i]+"_solid",
                               2.5*cm, 2.5*cm, 25.0*cm);               // half-lengths
        G4LogicalVolume* log = new G4LogicalVolume(box, matScintillator,
                               Scintillator_name_list[i]+name+"Logic");
        logicScintillators.push_back(log);
        physScint[i] = new G4PVPlacement(rot, translation, log,
                       Scintillator_name_list[i]+name, logicWorld, false, i, true);
    }

    // --------------------------------------------------------------
    // 4. Al foil – open at the SiPM face (z = -25 cm)
    // --------------------------------------------------------------
    const G4double foilThick = 0.0016*cm;               // 0.016 mm
    for (int i=0;i<nFoil;i++) {
        G4Box* outer = new G4Box("outerFoil",
                                 2.5*cm+foilThick, 2.5*cm+foilThick, 25.0*cm);
        G4Box* inner = new G4Box("innerFoil",
                                 2.5*cm, 2.5*cm, 25.0*cm-0.1*cm);   // shortened 1 mm
        G4ThreeVector shift(0,0,+0.1*cm);                  // open at -z
        G4SubtractionSolid* foil = new G4SubtractionSolid(
                Tapflon_name_list[i]+"_solid", outer, inner, nullptr, shift);

        G4LogicalVolume* log = new G4LogicalVolume(foil, matWrapping,
                               Tapflon_name_list[i]+name+"Logic");
        logicTapflon.push_back(log);
        physFoil[i] = new G4PVPlacement(rot, translation, log,
                       Tapflon_name_list[i]+name, logicWorld, false, i, true);
    }

    // --------------------------------------------------------------
    // 5. 4-sided Al protection (Bottom/Left/Right/Top)
    // --------------------------------------------------------------
    const G4double alThick = 1.0*cm;
    for (int i=0;i<nProt;i++) {
        G4Box* box = nullptr;
        G4ThreeVector localPos(0,0,0);

        if (i==0) {                                   // Bottom
            box = new G4Box(Protection_name_list[i]+"_solid",
                            2.5*cm+foilThick, alThick/2, 25.0*cm);
            localPos = G4ThreeVector(0,
                     -(2.5*cm+foilThick+alThick/2), 0);
        }
        else if (i==1) {                              // Left
            box = new G4Box(Protection_name_list[i]+"_solid",
                            alThick/2, 2.5*cm+foilThick, 25.0*cm);
            localPos = G4ThreeVector(-(2.5*cm+foilThick+alThick/2),0,0);
        }
        else if (i==2) {                              // Right
            box = new G4Box(Protection_name_list[i]+"_solid",
                            alThick/2, 2.5*cm+foilThick, 25.0*cm);
            localPos = G4ThreeVector(+(2.5*cm+foilThick+alThick/2),0,0);
        }
        else {                                        // Top
            box = new G4Box(Protection_name_list[i]+"_solid",
                            2.5*cm+foilThick, alThick/2, 25.0*cm);
            localPos = G4ThreeVector(0,
                     +(2.5*cm+foilThick+alThick/2), 0);
        }

        G4LogicalVolume* log = new G4LogicalVolume(box, matAl,
                               Protection_name_list[i]+name+"Logic");
        logicProtection.push_back(log);
        physProt[i] = new G4PVPlacement(rot,
                       translation + (*rot)(localPos),
                       log, Protection_name_list[i]+name,
                       logicWorld, false, i, true);
    }

    // --------------------------------------------------------------
    // 6. Acrylic plate – **same thickness as SiPM** (1.64 mm)
    // --------------------------------------------------------------
    const G4double plateThick = 0.164*cm;               // 1.64 mm
    for (int i=0;i<nAcry;i++) {
        // ---- base plate -------------------------------------------------
        G4Box* base = new G4Box(Acrylic_name_list[i]+"_base",
                                2.5*cm, 2.5*cm, plateThick/2);

        // ---- holes exactly the size of the SiPMs -----------------------
        G4Box* hole = new G4Box("hole",
                                0.209*cm, 0.209*cm, plateThick/2 + 0.01*mm); // tiny overlap for subtraction

        const G4double holeXY[4][2] = {
            {-1.25*cm, 1.25*cm},
            { 1.25*cm, 1.25*cm},
            {-1.25*cm,-1.25*cm},
            { 1.25*cm,-1.25*cm}
        };

        G4VSolid* acry = base;
        for (int j=0;j<4;j++) {
            G4ThreeVector hp(holeXY[j][0], holeXY[j][1], 0);
            acry = new G4SubtractionSolid(
                       Acrylic_name_list[i]+"_h"+std::to_string(j),
                       acry, hole, nullptr, hp);
        }

        G4LogicalVolume* log = new G4LogicalVolume(acry, matAcrylic,
                               Acrylic_name_list[i]+name+"Logic");
        logicAcrylic.push_back(log);

        // ---- place **exactly on the CsI face** (no air) ----------------
        G4ThreeVector posAcry(0,0,-25.0*cm - plateThick/2);
        physAcry[i] = new G4PVPlacement(rot,
                       translation + (*rot)(posAcry),
                       log, Acrylic_name_list[i]+name,
                       logicWorld, false, i, true);
    }

    // --------------------------------------------------------------
    // 7. 4 SiPMs – **embedded in the acrylic holes** (same thickness)
    // --------------------------------------------------------------
    for (int i=0;i<nSiPM;i++) {
        // half-length = 0.82 mm  → full thickness = 1.64 mm
        G4Box* box = new G4Box(SiPM_name_list[i]+"_solid",
                               0.209*cm, 0.209*cm, 0.082*cm);

        G4LogicalVolume* log = new G4LogicalVolume(box, matSiPM,
                               SiPM_name_list[i]+name+"Logic");
        logicCalorimeter = log;                 // for SD
        logicSiPM.push_back(log);

        const G4double holeXY[4][2] = {
            //{-1.25*cm, 1.25*cm},
            //{ 1.25*cm, 1.25*cm},
            {-1.25*cm,-1.25*cm},
            { 1.25*cm,-1.25*cm}
        };
        // centre of the hole = centre of the SiPM
        G4ThreeVector posSiPM(holeXY[i][0], holeXY[i][1],
                              -25.0*cm - plateThick/2);

        physSiPM[i] = new G4PVPlacement(rot,
                       translation + (*rot)(posSiPM),
                       log, SiPM_name_list[i]+name,
                       logicWorld, false, i, true);
    }

    // --------------------------------------------------------------
    // 8. Back Al Foil (behind Acrylic + SiPM layer)
    // --------------------------------------------------------------
    const G4double backFoilThick = 0.0016*cm;  // 0.016 mm
    G4Box* backFoilBox = new G4Box("BackFoil_solid",
                                   2.5*cm, 2.5*cm, backFoilThick/2);

    G4LogicalVolume* logicBackFoil = new G4LogicalVolume(backFoilBox, matWrapping,
                                                        "BackFoil" + name + "Logic");

    // Place it directly behind the acrylic layer (no gap)
    G4ThreeVector backFoilPos(0, 0, -25.0*cm - plateThick - backFoilThick/2);
    G4VPhysicalVolume* physBackFoil = new G4PVPlacement(
        rot, translation + (*rot)(backFoilPos),
        logicBackFoil, "BackFoil" + name, logicWorld, false, 0, true);

    // --------------------------------------------------------------
    // 9. Optical border surfaces (direction matters!)
    // --------------------------------------------------------------
    for (int i=0;i<nScint;i++) {
        // CsI → Acrylic (direct contact)
        new G4LogicalBorderSurface("CsI_Acrylic",
                                   physScint[i], physAcry[0], surfCsI_SiPM);

        // Acrylic → each SiPM (direct contact)
        for (int j=0;j<nSiPM;j++) {
            new G4LogicalBorderSurface("Acrylic_SiPM",
                                       physAcry[0], physSiPM[j], surfCsI_SiPM);
        }

        // CsI → Al foil (reflective)
        for (int j=0;j<nFoil;j++) {
            new G4LogicalBorderSurface("CsI_AlFoil",
                                       physScint[i], physFoil[j], surfCsI_AlFoil);
        }
        
        // Acrylic → Back Foil (reflects stray light back into SiPMs)
        new G4LogicalBorderSurface("Acrylic_BackFoil",
                                   physAcry[0], physBackFoil, surfCsI_AlFoil);
    }
}
void MyDetectorConstruction::ConstructCalorimeter_unit_3d(G4ThreeVector translation, G4double angle, G4String name){
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angle);
    rotation->rotateZ(30.0*deg); //Remove this line if no self rotation
    std::string Scintillator_name_list[] = {"Hexagonal/UntitledPrism12.2mm"};
    std::string SiPM_name_list[] = {"Hexagonal/UntitledSiPM1_12.2mm", "Hexagonal/UntitledSiPM2_12.2mm",
                                    "Hexagonal/UntitledSiPM3_12.2mm", "Hexagonal/UntitledSiPM4_12.2mm"};
    std::string Tapflon_name_list[] = {"Hexagonal/UntitledTape12.2mm"};
    int Size_of_Scintillator_name_list = sizeof(Scintillator_name_list)/sizeof(std::string);
    int Size_of_SiPM_name_list = sizeof(SiPM_name_list)/sizeof(std::string);
    int Size_of_Tapflon_name_list = sizeof(Tapflon_name_list)/sizeof(std::string);

    std::vector<G4VPhysicalVolume*> physScintillators(Size_of_Scintillator_name_list);
    std::vector<G4VPhysicalVolume*> physTapflon(Size_of_Tapflon_name_list);
    std::vector<G4VPhysicalVolume*> physSiPM(Size_of_SiPM_name_list);
    for (int i = 0; i < Size_of_Scintillator_name_list; i++) {
        std::string name_scint = Scintillator_name_list[i];
        auto scintillator = CADMesh::TessellatedMesh::FromSTL(name_scint + ".stl");
        scintillator->SetScale(1.0);
        auto Scintillator = scintillator->GetSolid();
        G4LogicalVolume* logicScintillator_pre = new G4LogicalVolume(Scintillator, matScintillator, name_scint+name + "Logic");
        logicScintillators.push_back(logicScintillator_pre);
        physScintillators[i] = new G4PVPlacement(rotation, translation, logicScintillator_pre, name_scint+name, logicWorld, false, i, true);    

        std::string name_SiPM = SiPM_name_list[i];
        auto scintillatorDet = CADMesh::TessellatedMesh::FromSTL(name_SiPM + ".stl");
        scintillatorDet->SetScale(1.0);
        auto ScintillatorDet = scintillatorDet->GetSolid();
        G4LogicalVolume* logicSiPM_pre = new G4LogicalVolume(ScintillatorDet, matSiPM, name_SiPM+name + "Logic");
		logicCalorimeter=logicSiPM_pre;
        logicSiPM.push_back(logicCalorimeter);
        physSiPM[i] = new G4PVPlacement(rotation, translation, logicCalorimeter, name_SiPM+name, logicWorld, false, i, true);    

        std::string name_Wrapping = Tapflon_name_list[i];
        auto scintillatorwrapping = CADMesh::TessellatedMesh::FromSTL(name_Wrapping + ".stl");
        scintillatorwrapping->SetScale(1.0);
        auto Scintillatorwrapping = scintillatorwrapping->GetSolid();
        G4LogicalVolume* logicTapflon_pre = new G4LogicalVolume(Scintillatorwrapping, matWrapping, name_Wrapping+name + "Logic");
        logicTapflon.push_back(logicTapflon_pre);
        physTapflon[i] = new G4PVPlacement(rotation, translation, logicTapflon_pre, name_Wrapping+name, logicWorld, false, i, true);    
    }
    for (int i=0; i<Size_of_Scintillator_name_list; i++) {
        new G4LogicalBorderSurface("CsI_SiPM_Border", physScintillators[i], physSiPM[i], surfCsI_SiPM);
        new G4LogicalBorderSurface("CsI_Teflon_Border", physScintillators[i], physTapflon[i], surfCsI_Teflon);
        //new G4LogicalBorderSurface("CsI_SiPM_Border_Reverse", physSiPM[i], physScintillators[i], surfCsI_SiPM);
        //new G4LogicalBorderSurface("CsI_Teflon_Border_Reverse", physTapflon[i], physScintillators[i], surfCsI_Teflon);
    }
}



void MyDetectorConstruction::ConstructCalorimeter() {
    // Place a single unit at origin
    if(is3DCalorimeter){
        // Generate cubic
        //int range=0;
        //G4double dist=0*mm;
        //int counter=0;
        //for(int j=0;j<=range;j++){
        //    for(int i=-range;i<=range;i++){
        //        for(int k=-range;k<=range;k++){
        //            G4String name_=to_string(i)+"_"+to_string(j)+"_"+to_string(k);
        //            G4double angle = 0 * deg;
        //            if(i==0&&j==0&&k==0){
        //                G4ThreeVector translation(0.*mm+(i*6.05*2)*mm, 0.*mm+(j*6.05*2)*mm, 0.*mm+(k*6.05*2)*mm);
        //                ConstructCalorimeter_unit(translation,angle,name_);
        //                counter+=1;
        //            }
        //            else{
        //                G4ThreeVector translation(0.*mm+(i*(6.05)*2+std::copysign(1.0f,i)*dist)*mm, 0.*mm+(j*(6.05)*2+std::copysign(1.0f,j)*dist)*mm, 0.*mm+(k*(6.05)*2+std::copysign(1.0f,k)*dist)*mm);
        //                ConstructCalorimeter_unit(translation,angle,name_);
        //            }       
        //
        //        }
        //    }
        //}
        // Generate hcc
        G4double apothem = (12.2+0.3)/2*std::sqrt(3.0)/2.0;  // Apothem (distance from center to flat side)
        G4double side_length = 2.0 * apothem;  // Side length
        G4double a1_x = side_length;  // Primitive vector 1 x-component
        G4double a1_y = 0.0;  // Primitive vector 1 y-component
        G4double a2_x = side_length / 2.0;  // Primitive vector 2 x-component
        G4double a2_y = side_length * std::sqrt(3.0) / 2.0;  // Primitive vector 2 y-component (sin(60°))

        int min_N = 13;  // Start from ring 1 for placing source
        int max_N = 13+(3-1);  // End at ring 6
        int count = 0;  // For unique naming

        for (int n1 = -max_N; n1 <= max_N; ++n1) {
            for (int n2 = std::max(-max_N, -n1 - max_N); n2 <= std::min(max_N, -n1 + max_N); ++n2) {
            // Calculate the "ring" distance from origin using the max norm
            int ring = std::max({std::abs(n1), std::abs(n2), std::abs(n1 + n2)});
            if (ring >= min_N && ring <= max_N) {
                // Calculate position using primitive vectors
                G4double x = n1 * a1_x + n2 * a2_x;
                G4double y = n1 * a1_y + n2 * a2_y;
                G4double z = 0.0;  // Adjust if prisms are offset along z

                G4ThreeVector translation(x, y, z);  // Units: assume bare numbers match your radius units
                G4double angle = 0.0*deg;  // No rotation; adjust if needed to align with prism definition
                G4String name = "calor_unit_" + std::to_string(count++);

                // Call your function to place the unit
                ConstructCalorimeter_unit_3d(translation, angle, name);
            }
            }
        }
    }
    else{
        ConstructCalorimeter_unit(G4ThreeVector(0,-(25*cm-(4*bare_source_radius)),(-1*cm+3.5*cm+disk_height_half+ring_height_half)), 90*deg, "");
    }
}
//Construct source
void MyDetectorConstruction::ConstructSource(){
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

void MyDetectorConstruction::ConstructLiquidScintillator(){
	G4Tubs* solidContainer = new G4Tubs("LiquidContainer", 0.*nm, container_radius, container_height_half, 0.*deg, 360.*deg);
	G4Tubs* solidLiquid = new G4Tubs("Liquid", 0.*nm, container_radius-container_thickness, container_height_half-container_thickness, 0.*deg, 360.*deg);

	logiContainer_F = new G4LogicalVolume(solidContainer, matContainer, "logiContainer");
	logiContainer_B = new G4LogicalVolume(solidContainer, matContainer, "logiContainer");
	logicLiquid_F = new G4LogicalVolume(solidLiquid, matLiquid, "logicLiquid");
	logicLiquid_B = new G4LogicalVolume(solidLiquid, matLiquid, "logicLiquid");

	G4double container_z_shift = ring_height_half+container_height_half;

	G4Translate3D transZ(G4ThreeVector(0.*m, 0.*m, container_z_shift));
	G4Rotate3D rotY_90_2(90.*2*deg, G4ThreeVector(0.,1.,0.));

	physContainer_F = new G4PVPlacement(transZ, logiContainer_F, "Container_F", logicWorld, false, 0, true);
	physContainer_B = new G4PVPlacement(rotY_90_2*transZ, logiContainer_B, "Container_B", logicWorld, false, 0, true);
	physLiquid_F = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, 0.*m), logicLiquid_F, "Liquid", logiContainer_F, false, 0, true);
	physLiquid_B = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, 0.*m), logicLiquid_B, "Liquid", logiContainer_B, false, 0, true);
}
// End Construct source