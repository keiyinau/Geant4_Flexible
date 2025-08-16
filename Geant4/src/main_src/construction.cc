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
	matWorld = Air; //Vacuum;

	// Set the default of each logical volume to be NULL so the sensitive detector selector can work well
	logicDetector_Shell = NULL;
	logicTPC = NULL;
	logicCalorimeter = NULL;

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
    G4MaterialPropertiesTable* mptAir = new G4MaterialPropertiesTable();
    mptAir->AddProperty("RINDEX", "Air");
    Air->SetMaterialPropertiesTable(mptAir);

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

	//NaCl
	matNaCl = new G4Material("NaCl", 2.16*g/cm3, 2);
	matNaCl->AddElement(nist->FindOrBuildElement("Na"), 1);
	matNaCl->AddElement(nist->FindOrBuildElement("Cl"), 4);
	// End NaCl
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
	matCsI->SetMaterialPropertiesTable(mptCsI);
	// End CsI

	// Define Tapflon(teflon) for wrapping
	matTeflon = nist->FindOrBuildMaterial("G4_TEFLON");
	std::vector<G4double> tapflon_reflectance_Energy, tapflon_reflectance_fractions;
	readAndProcessData_Energy_txt("teflon_Reflectance-modified.txt", tapflon_reflectance_Energy, tapflon_reflectance_fractions);
    std::vector<G4double> tapflon_refraction_Energy, tapflon_refraction_Index;
	readAndProcessData_Energy_txt("Refraction_Index_Teflon_Gray.txt", tapflon_refraction_Energy, tapflon_refraction_Index);
	G4MaterialPropertiesTable* mptTeflon = new G4MaterialPropertiesTable();
	mptTeflon->AddProperty("REFLECTIVITY", tapflon_reflectance_Energy, tapflon_reflectance_fractions,tapflon_reflectance_fractions.size());
    mptTeflon->AddProperty("RINDEX", tapflon_refraction_Energy, tapflon_refraction_Index,tapflon_refraction_Index.size());

    matTeflon->SetMaterialPropertiesTable(mptTeflon);
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
	matSi->SetMaterialPropertiesTable(mptSi);

    // CsI-SiPM (dielectric-dielectric interface)
    surfCsI_SiPM = new G4OpticalSurface("CsI_SiPM_Surface");
    surfCsI_SiPM->SetType(dielectric_dielectric);
    surfCsI_SiPM->SetModel(unified); // Glisur for smooth dielectric interface
    surfCsI_SiPM->SetFinish(polished);
	// End SiPM

    // CsI-Teflon (reflective surface)
    surfCsI_Teflon = new G4OpticalSurface("CsI_Teflon_Surface");
    surfCsI_Teflon->SetType(dielectric_metal); // Teflon as reflective surface
    surfCsI_Teflon->SetModel(unified);
    surfCsI_Teflon->SetFinish(polished);
    //End surface



    std::cout<<"==========================="<<std::endl;
    std::cout<<"Printing the material properties of CsI"<<std::endl;
    mptCsI->DumpTable();
    std::cout<<"==========================="<<std::endl;
    std::cout<<"==========================="<<std::endl;
    std::cout<<"Printing the material properties of Teflon"<<std::endl;
    mptTeflon->DumpTable();
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
	G4double xWorld = 1*m;
	G4double yWorld = 1*m;
	G4double zWorld = 1*m;

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
    std::cout<<"==========================="<<std::endl;
    std::cout<<"Test if there are overlap"<<std::endl;
    physWorld->CheckOverlaps();
    std::cout<<"==========================="<<std::endl;

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
        for(int i=0; i < logicSiPM.size(); i++) {
            logicSiPM[i]->SetSensitiveDetector(calorimeter);
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
	G4double shell_thickness = 3.*cm;//1.*nm;
	G4double inner_radius =0.*cm;// 25.*cm+80.*cm;
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
void MyDetectorConstruction::ConstructCalorimeter_unit(G4ThreeVector translation, G4double angle, G4String name){
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angle);

    std::string Scintillator_name_list[] = {"Scintillator"};
    std::string SiPM_name_list[] = {"SiPM"};
    std::string Tapflon_name_list[] = {"Wrapping"};
    int Size_of_Scintillator_name_list = sizeof(Scintillator_name_list)/sizeof(std::string);
    int Size_of_SiPM_name_list = sizeof(SiPM_name_list)/sizeof(std::string);
    int Size_of_Tapflon_name_list = sizeof(Tapflon_name_list)/sizeof(std::string);

    std::vector<G4VPhysicalVolume*> physScintillators(Size_of_Scintillator_name_list);
    std::vector<G4VPhysicalVolume*> physTapflon(Size_of_Tapflon_name_list);
    std::vector<G4VPhysicalVolume*> physSiPM(Size_of_SiPM_name_list);
    for (int i = 0; i < Size_of_Scintillator_name_list; i++) {
        std::string name_scint = Scintillator_name_list[i];
        auto scintillator = CADMesh::TessellatedMesh::FromSTL(name_scint + ".stl");
        scintillator->SetScale(10.0);
        auto Scintillator = scintillator->GetSolid();
        G4LogicalVolume* logicScintillator_pre = new G4LogicalVolume(Scintillator, matCsI, name_scint+name + "Logic");
        logicScintillators.push_back(logicScintillator_pre);
        physScintillators[i] = new G4PVPlacement(rotation, translation, logicScintillator_pre, name_scint+name, logicWorld, false, i, true);    

        std::string name_SiPM = SiPM_name_list[i];
        auto scintillatorDet = CADMesh::TessellatedMesh::FromSTL(name_SiPM + ".stl");
        scintillatorDet->SetScale(10.0);
        auto ScintillatorDet = scintillatorDet->GetSolid();
        G4LogicalVolume* logicSiPM_pre = new G4LogicalVolume(ScintillatorDet, matSi, name_SiPM+name + "Logic");
		logicCalorimeter=logicSiPM_pre;
        logicSiPM.push_back(logicCalorimeter);
        physSiPM[i] = new G4PVPlacement(rotation, translation, logicCalorimeter, name_SiPM+name, logicWorld, false, i, true);    

        std::string name_Wrapping = Tapflon_name_list[i];
        auto scintillatorwrapping = CADMesh::TessellatedMesh::FromSTL(name_Wrapping + ".stl");
        scintillatorwrapping->SetScale(10.0);
        auto Scintillatorwrapping = scintillatorwrapping->GetSolid();
        G4LogicalVolume* logicTapflon_pre = new G4LogicalVolume(Scintillatorwrapping, matTeflon, name_Wrapping+name + "Logic");
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
    int range=30;
    G4double dist=0*mm;
    int counter=0;
    for(int j=0;j<=range;j++){
        for(int i=-range;i<=range;i++){
            for(int k=-range;k<=range;k++){
                G4String name_=to_string(i)+"_"+to_string(j)+"_"+to_string(k);
                G4double angle = 90 * deg;
                if(i==0&&j==0&&k==0){
                    G4ThreeVector translation(0.*mm+(i*6.05*2)*mm, 0.*mm+(j*6.05*2)*mm, 0.*mm+(k*6.05*2)*mm);
                    ConstructCalorimeter_unit(translation,angle,name_);
                    counter+=1;
                }
                else{
                    G4ThreeVector translation(0.*mm+(i*(6.05)*2+std::copysign(1.0f,i)*dist)*mm, 0.*mm+(j*(6.05)*2+std::copysign(1.0f,j)*dist)*mm, 0.*mm+(k*(6.05)*2+std::copysign(1.0f,k)*dist)*mm);
                    ConstructCalorimeter_unit(translation,angle,name_);
                }       

            }
        }

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