#include "construction.hh"
#include "CADMesh.hh"
MyDetectorConstruction::MyDetectorConstruction() {
	// Define required materials
    logicOptical=true;
    fLightYield=1000;		// Light yield of scintillator, in photons/MeV. This is a made-up value for demonstration; adjust based on actual material properties.

	DefineMaterials();


    coordinate_name="coordinates.txt";
	isDetector_Shell = false;
	isSource=false;
	isTPC = false;
	isCalorimeter = true;
    isLiquid=false;
    is3DCalorimeter=true;
	// Set the material for each logical volume
	matWorld = Air; //Vacuum;
    matLiquid=matAcrylic;
    matContainer=matAcrylic;
    matScintillator=matWater;
    matSiPM=matSi;
    matWrapping=matAl;
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
    container_radius = 1*cm;				//All variables of container are made up.
	container_height_half = 0.5/2*cm;
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
bool MyDetectorConstruction::readAndProcessData_Nonproportionality(const std::string& filename, 
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
        emission_Energy.push_back(p.first*1e-3);
        emission_fractions.push_back(p.second);
    }

    return true;
}
G4String MyDetectorConstruction::file_name = "";

void MyDetectorConstruction::DefineMaterials() {
    G4NistManager* nist = G4NistManager::Instance();

    // ==========================================
    // 1. World Material (Air) - Air Coupling 必需
    // ==========================================
    Air = nist->FindOrBuildMaterial("G4_AIR");
    G4MaterialPropertiesTable* mptAir = new G4MaterialPropertiesTable();
    
    std::vector<G4double> Air_absorption_Energy, Air_absorption_Index;
    readAndProcessData_Energy_cm_txt("AbsorptionLength_Air.txt", Air_absorption_Energy, Air_absorption_Index);
    
    // 顯式定義 Air 的 RINDEX 為 1.0 (極度重要，確保光子穿透表面後不被刪除)
    G4double airEnergies[] = { 2.034*eV, 2.384*eV, 2.755*eV, 3.100*eV }; 
    G4double rIndexAir[] = { 1.0, 1.0, 1.0, 1.0 }; 
    mptAir->AddProperty("RINDEX", airEnergies, rIndexAir, 4);
    mptAir->AddProperty("ABSLENGTH", Air_absorption_Energy, Air_absorption_Index, Air_absorption_Index.size());

    // ==========================================
    // 2. Radioactive Source Materials
    // ==========================================
    // NaCl (Active Source) - 修正比例為 1:1
    matNaCl = new G4Material("NaCl", 2.16*g/cm3, 2);
    matNaCl->AddElement(nist->FindOrBuildElement("Na"), 1);
    matNaCl->AddElement(nist->FindOrBuildElement("Cl"), 1);

    // Ti Foil
    G4Element* Ti = nist->FindOrBuildElement("Ti");
    matTi = new G4Material("Ti_", 4.507*g/cm3, 1);
    matTi->AddElement(Ti, 1.);

    // ==========================================
    // 3. Setup Components & Elements
    // ==========================================
    G4Element* elH  = nist->FindOrBuildElement("H");
    G4Element* elC  = nist->FindOrBuildElement("C");
    G4Element* elO  = nist->FindOrBuildElement("O");
    G4Element* elSi = nist->FindOrBuildElement("Si");
    G4Element* elCl = nist->FindOrBuildElement("Cl");

    G4double photonEnergy[] = { 2.0*eV, 2.5*eV, 3.0*eV, 3.5*eV };
    const G4int nEntries = sizeof(photonEnergy) / sizeof(G4double);

    // Cuvette (Quartz / Fused Silica)
    matCuvette = new G4Material("matCuvette", 2.20*g/cm3, 2);
    matCuvette->AddElement(elSi, 1);
    matCuvette->AddElement(elO, 2);
    G4MaterialPropertiesTable* mptCuvette = new G4MaterialPropertiesTable();
    G4double rindexCuvette[] = { 1.46, 1.46, 1.46, 1.46 }; 
    mptCuvette->AddProperty("RINDEX", photonEnergy, rindexCuvette, nEntries);
    G4double absCuvette[] = { 10.0*m, 10.0*m, 10.0*m, 10.0*m }; // 假設為高透光，吸收長度 10 米
    mptCuvette->AddProperty("ABSLENGTH", photonEnergy, absCuvette, nEntries);

    // matPVC (Transparent Polyvinyl Chloride)
    // ==========================================
    matPVC = new G4Material("matPVC", 1.35*g/cm3, 3);
    matPVC->AddElement(elC, 2);
    matPVC->AddElement(elH, 3);
    matPVC->AddElement(elCl, 1);

    G4MaterialPropertiesTable* mptPVC = new G4MaterialPropertiesTable();

    // 折射率 (RINDEX): 透明 PVC 一般約為 1.52 - 1.54
    G4double rindexPVC[] = { 1.53, 1.53, 1.53, 1.53 };
    mptPVC->AddProperty("RINDEX", photonEnergy, rindexPVC, nEntries);

    // 吸收長度 (ABSLENGTH): 模擬透明 PVC 的有限透光度，設為 2.0 米
    G4double absPVC[] = { 2.0*m, 2.0*m, 2.0*m, 2.0*m };
    mptPVC->AddProperty("ABSLENGTH", photonEnergy, absPVC, nEntries);

    matPVC->SetMaterialPropertiesTable(mptPVC);

    // matPhotopolymer (3D Printing Resin - Clear Acrylic Base)
    // ==========================================
    matPhotopolymer = new G4Material("matPhotopolymer", 1.20*g/cm3, 3);
    matPhotopolymer->AddElement(elC, 5);
    matPhotopolymer->AddElement(elH, 8);
    matPhotopolymer->AddElement(elO, 2);

    G4MaterialPropertiesTable* mptPhotopolymer = new G4MaterialPropertiesTable();

    // 折射率 (RINDEX): 標準光固化樹脂 (PMMA base) 約為 1.49
    G4double rindexPhotopolymer[] = { 1.49, 1.49, 1.49, 1.49 };
    mptPhotopolymer->AddProperty("RINDEX", photonEnergy, rindexPhotopolymer, nEntries);

    // 吸收長度 (ABSLENGTH): 透明 3D 打印樹脂透光度有限，設為 3.0 cm
    // (如果你用的是灰色/黑色不透明樹脂，請不要為此物料設定 RINDEX 和 ABSLENGTH)
    G4double absPhotopolymer[] = { 3.0*cm, 3.0*cm, 3.0*cm, 3.0*cm };
    mptPhotopolymer->AddProperty("ABSLENGTH", photonEnergy, absPhotopolymer, nEntries);

    matPhotopolymer->SetMaterialPropertiesTable(mptPhotopolymer);

    // ==========================================
    // 4. Scintillators
    // ==========================================
    
    // --- Liquid Scintillator (Generator) ---
    matGenerator = new G4Material("matGenerator", 0.867*g/cm3, 2);
    matGenerator->AddElement(elC, 7);
    matGenerator->AddElement(elH, 8);
    G4MaterialPropertiesTable* mptGenerator = new G4MaterialPropertiesTable();
    
    G4double rindexGenerator[] = { 1.50, 1.50, 1.50, 1.50 }; 
    mptGenerator->AddProperty("RINDEX", photonEnergy, rindexGenerator, nEntries);
    
    // 載入 Daya Bay LAB 發光數據
    //std::vector<G4double> LS_emission_Energy, LS_emission_fractions;
    //readAndProcessData_Energy("LAB_DayaBay_Normalized.csv", LS_emission_Energy, LS_emission_fractions);
    //mptGenerator->AddProperty("SCINTILLATIONCOMPONENT1", LS_emission_Energy, LS_emission_fractions, LS_emission_fractions.size());
    //mptGenerator->AddConstProperty("SCINTILLATIONYIELD", fLightYield/MeV);
    //mptGenerator->AddConstProperty("SCINTILLATIONYIELD1", 1.0); 
    //mptGenerator->AddConstProperty("RESOLUTIONSCALE", 1.0);
    //mptGenerator->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 5.0 * ns);
    // 你可以根據 Daya Bay 或相關文獻調整此數值，這裡暫定為典型的 5 米
    G4double absGenerator[] = { 5.0*m, 5.0*m, 5.0*m, 5.0*m }; 
    mptGenerator->AddProperty("ABSLENGTH", photonEnergy, absGenerator, nEntries);
    // ==========================================
    // --- 4. Plastic Scintillator (PVT - EJ-200 / BC-408 Equivalent) ---
    // ==========================================
    matPlasticScintillator = new G4Material("matPlasticScintillator", 1.023*g/cm3, 2);
    matPlasticScintillator->AddElement(elC, 9);
    matPlasticScintillator->AddElement(elH, 10);
    
    G4MaterialPropertiesTable* mptPlasticScint = new G4MaterialPropertiesTable();
    
    // 1. 定義能量陣列 (必須由低能量/長波長 排列至 高能量/短波長)
    // 對應波長: 500nm, 475nm, 450nm, 425nm(Peak), 410nm, 400nm, 380nm
    G4double pvtEnergies[] = { 2.48*eV, 2.61*eV, 2.75*eV, 2.92*eV, 3.02*eV, 3.10*eV, 3.26*eV };
    const G4int numPvtEntries = sizeof(pvtEnergies)/sizeof(G4double);
    
    // 2. 基礎光學屬性: 折射率 (RINDEX) 與 吸收長度 (ABSLENGTH)
    // EJ-200 嘅 Refractive Index 係 1.58，Bulk Attenuation Length 係 380 cm
    G4double rindexPvt[] = { 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58 }; 
    G4double absPvt[]    = { 380.0*cm, 380.0*cm, 380.0*cm, 380.0*cm, 380.0*cm, 380.0*cm, 380.0*cm };
    
    mptPlasticScint->AddProperty("RINDEX", pvtEnergies, rindexPvt, numPvtEntries);
    mptPlasticScint->AddProperty("ABSLENGTH", pvtEnergies, absPvt, numPvtEntries);
    
    // 3. 發光頻譜 (Scintillation Emission Spectrum)
    // 根據 EJ-200 規格曲線設定嘅相對強度 (Peak 位於 425nm / 2.92 eV)
    G4double emissionPvt[] = { 0.05, 0.20, 0.50, 1.00, 0.60, 0.10, 0.00 };
    mptPlasticScint->AddProperty("SCINTILLATIONCOMPONENT1", pvtEnergies, emissionPvt, numPvtEntries);
    
    // 4. 發光產量與時間常數 (Yield & Time Constants)
    mptPlasticScint->AddConstProperty("SCINTILLATIONYIELD", 10000.0/MeV); // 每 MeV 產生 10,000 光子
    mptPlasticScint->AddConstProperty("SCINTILLATIONYIELD1", 1.0);        // 100% 權重分配比呢個 Component
    mptPlasticScint->AddConstProperty("RESOLUTIONSCALE", 1.0);            // 預設統計分佈
    mptPlasticScint->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.1 * ns); // 衰減時間 (Decay Time)
    mptPlasticScint->AddConstProperty("SCINTILLATIONRISETIME1", 0.9 * ns);     // 上升時間 (Rise Time)
    
    // 將 Properties Table 綁定至 Material
    matPlasticScintillator->SetMaterialPropertiesTable(mptPlasticScint);
    // ==========================================
    // 5. Sensors (SiPM)
    // ==========================================
    matSi = nist->FindOrBuildMaterial("G4_Si");
    std::vector<G4double> Si_reflectance_Energy, Si_reflectance_fractions;
    readAndProcessData_Energy_txt("Reflectance_Si.txt", Si_reflectance_Energy, Si_reflectance_fractions);
    std::vector<G4double> Si_transmission_Energy, Si_rtransmission_Index;
    readAndProcessData_txt("transmittance_Si.txt", Si_transmission_Energy, Si_rtransmission_Index);
    std::vector<G4double> Si_refraction_Energy, Si_refraction_Index;
    readAndProcessData_txt("RefractiveIndexINFO_Si.txt", Si_refraction_Energy, Si_refraction_Index);

    G4MaterialPropertiesTable* mptSi = new G4MaterialPropertiesTable();
    mptSi->AddProperty("REFLECTIVITY", Si_reflectance_Energy, Si_reflectance_fractions, Si_reflectance_fractions.size());
    mptSi->AddProperty("TRANSMITTANCE", Si_transmission_Energy, Si_rtransmission_Index, Si_rtransmission_Index.size());  
    mptSi->AddProperty("RINDEX", Si_refraction_Energy, Si_refraction_Index, Si_refraction_Energy.size());    

    // ==========================================
    // 6. Optical Physics Toggle
    // ==========================================`
    // ==========================================
    if(logicOptical){
        Air->SetMaterialPropertiesTable(mptAir);
        matCuvette->SetMaterialPropertiesTable(mptCuvette);
        matGenerator->SetMaterialPropertiesTable(mptGenerator);
        matPlasticScintillator->SetMaterialPropertiesTable(mptPlasticScint);
        matSi->SetMaterialPropertiesTable(mptSi);
        matPhotopolymer->SetMaterialPropertiesTable(mptPhotopolymer);
        
        // [新增] 載入透明 PVC 嘅光學屬性
        matPVC->SetMaterialPropertiesTable(mptPVC); 
    }

    // ==========================================
    // 7. Debug Dumps
    // ==========================================
    std::cout<<"==========================="<<std::endl;
    std::cout<<"Printing the material properties of Liquid Scintillator (Generator)"<<std::endl;
    mptGenerator->DumpTable();
    std::cout<<"==========================="<<std::endl;
    std::cout<<"Printing the material properties of Plastic Scintillator"<<std::endl;
    mptPlasticScint->DumpTable();
    std::cout<<"==========================="<<std::endl;
    std::cout<<"Printing the material properties of Air"<<std::endl;
    mptAir->DumpTable();
    std::cout<<"==========================="<<std::endl;
    std::cout<<"Printing the material properties of Si (Sensors)"<<std::endl;
    mptSi->DumpTable();
    std::cout<<"==========================="<<std::endl;
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
    fMessenger->DeclareProperty("setDetectorCoordinate",coordinate_name,"Set the coordinate file that is using");
    fMessenger->DeclareMethod("setLightYield", &MyDetectorConstruction::SetLightYield, "Set Light Yield in photons/MeV");
}
// Construct All physical volumes
G4VPhysicalVolume* MyDetectorConstruction::Construct() {
	G4double xWorld = 0.1*m;
	G4double yWorld = 0.1*m;
	G4double zWorld = 0.1*m;

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
		logicDetector_Shell->SetSensitiveDetector(detect_edep);
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
	G4double shell_thickness = 2.*cm;//1.*nm;
	G4double inner_radius =5.0*cm;// 25.*cm+80.*cm;
	G4double outer_radius = inner_radius + shell_thickness;
	G4Sphere* solidDetector_Shell = new G4Sphere("solidDetector_Shell", inner_radius, outer_radius, 0.*deg, 360.*deg, 0.*deg, 360.*deg);
	logicDetector_Shell = new G4LogicalVolume(solidDetector_Shell, matLiquid, "logicDetector_Shell");
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
    // 2. LYSO crystal 2×2×20 mm (half-sizes: 1×1×10 mm)
    // --------------------------------------------------------------
    G4Box* crystalSolid = new G4Box("LYSO_solid", 1.*mm, 1.*mm, 10.*mm);
    G4LogicalVolume* logicCrystal = new G4LogicalVolume(crystalSolid, matLYSO,
                                                        "LYSO" + name + "Logic");
    logicScintillators.push_back(logicCrystal);

    G4VPhysicalVolume* physCrystal = new G4PVPlacement(rot, translation, logicCrystal,
                                                       "LYSO" + name, logicWorld, false, 0, true);

    // --------------------------------------------------------------
    // 3. Teflon wrapping (5 faces, open at SiPM end)
    // --------------------------------------------------------------
    const G4double foilThick = 0.1*mm;
    G4Box* outer = new G4Box("outerFoil", 1.*mm + foilThick, 1.*mm + foilThick, 10.*mm);
    G4Box* inner = new G4Box("innerFoil", 1.*mm, 1.*mm, 10.*mm - 0.1*mm);
    G4ThreeVector shift(0, 0, +0.1*mm);   // open at -z (SiPM side)
    G4SubtractionSolid* foilSolid = new G4SubtractionSolid("Teflon_solid", outer, inner, nullptr, shift);

    G4LogicalVolume* logicTeflon = new G4LogicalVolume(foilSolid, matWrapping,
                                                       "Teflon" + name + "Logic");
    logicTapflon.push_back(logicTeflon);

    G4VPhysicalVolume* physTeflon = new G4PVPlacement(rot, translation, logicTeflon,
                                                      "Teflon" + name, logicWorld, false, 0, true);

    // --------------------------------------------------------------
    // 4. Optical grease layer (200 µm)
    // --------------------------------------------------------------
    G4double greaseThick = 0.2*mm;
    G4Box* greaseSolid = new G4Box("Grease_solid", 1.*mm, 1.*mm, greaseThick/2);
    G4Material* greaseMat = new G4Material("OpticalGrease", 1.05*g/cm3, 1);
    greaseMat->AddElement(G4Element::GetElement("C"), 0.6);
    G4LogicalVolume* logicGrease = new G4LogicalVolume(greaseSolid, greaseMat, "Grease" + name + "Logic");

    G4ThreeVector greasePos(0, 0, -10.*mm - greaseThick/2);
    G4VPhysicalVolume* physGrease = new G4PVPlacement(rot, translation + (*rot)(greasePos),
                                                      logicGrease, "Grease" + name, logicWorld, false, 0, true);

    // --------------------------------------------------------------
    // 5. One SiPM attached at the grease end
    // --------------------------------------------------------------
    G4Box* sipmSolid = new G4Box("SiPM_solid", 1.*mm, 1.*mm, 0.082*cm);  // 1.64 mm thick
    G4LogicalVolume* logicSiPMs = new G4LogicalVolume(sipmSolid, matSiPM, "SiPM" + name + "Logic");
    logicCalorimeter = logicSiPMs;   // for your SensitiveDetector

    G4ThreeVector sipmPos(0, 0, -10.*mm - greaseThick - 0.082*cm);
    G4VPhysicalVolume* physSiPM = new G4PVPlacement(rot, translation + (*rot)(sipmPos),
                                                    logicSiPMs, "SiPM" + name, logicWorld, false, 0, true);
    logicSiPM.push_back(logicCalorimeter);
    // --------------------------------------------------------------
    // 6. Optical border surfaces (direction matters!)
    // --------------------------------------------------------------
    // Crystal → Grease
    new G4LogicalBorderSurface("Crystal-Grease", physCrystal, physGrease, surfCrystalGrease);

    // Grease → SiPM
    new G4LogicalBorderSurface("Grease-SiPM", physGrease, physSiPM, surfCsI_SiPM);

    // Crystal → Teflon (reflective)
    new G4LogicalBorderSurface("Crystal-Teflon", physCrystal, physTeflon, surfCsI_Teflon);
}
void MyDetectorConstruction::ConstructCalorimeter_unit_3d(G4ThreeVector translation, G4String name, G4double rotateX, G4double rotateY, G4double rotateZ){
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(rotateX);
    rotation->rotateY(rotateY);
    rotation->rotateZ(rotateZ); 
    

    // === 1. DEFINE STL FILE PREFIX ===
    std::string prefix = "Generator/";
    
    // === 2.Radioactive Source ===
    auto activesourceMesh = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_Na22Source_1_Na22Active.stl");
    activesourceMesh->SetScale(1.0);
    G4LogicalVolume* logicactivesourceMesh = new G4LogicalVolume(activesourceMesh->GetSolid(), matNaCl, "logicactivesourceNa22" + name + "_Logic");
    G4VPhysicalVolume* physactivesource = new G4PVPlacement(rotation, translation, logicactivesourceMesh, "activesourceNa22" + name, logicWorld, false, 0, true);
   
    auto FoilsourceMesh = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_Na22Source_1_TiFoil.stl");
    FoilsourceMesh->SetScale(1.0);
    G4LogicalVolume* logicFoilsourceMesh = new G4LogicalVolume(FoilsourceMesh->GetSolid(), matTi, "logicFoilsourceTiRing" + name + "_Logic");
    G4VPhysicalVolume* physFoilsource = new G4PVPlacement(rotation, translation, logicFoilsourceMesh, "FoilsourceTiRing" + name, logicWorld, false, 0, true);
   
    auto PVCsourceMesh = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_SourcePVC_1_Body1.stl");
    PVCsourceMesh->SetScale(1.0);
    G4LogicalVolume* logicPVCsourceMesh = new G4LogicalVolume(PVCsourceMesh->GetSolid(), matPVC, "logicPVCsourcePVC" + name + "_Logic");
    G4VPhysicalVolume* physPVCsource = new G4PVPlacement(rotation, translation, logicPVCsourceMesh, "PVCsourcePVC" + name, logicWorld, false, 0, true);
   
    // === Plastic scintillator ===
    auto plasticscint1Mesh = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_PlasticScint_1_Body1.stl");
    plasticscint1Mesh->SetScale(1.0);
    G4LogicalVolume* logicPlasticScint1 = new G4LogicalVolume(plasticscint1Mesh->GetSolid(), matPlasticScintillator, "PlasticScint_" + name + "_Logic");
    G4VPhysicalVolume* physPlasticScint1 = new G4PVPlacement(rotation, translation, logicPlasticScint1, "PlasticScint_" + name, logicWorld, false, 0, true);
    logicScintillators.push_back(logicPlasticScint1);

    auto plasticscint2Mesh = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_PlasticScint_2_Body1.stl");
    plasticscint2Mesh->SetScale(1.0);
    G4LogicalVolume* logicPlasticScint2 = new G4LogicalVolume(plasticscint2Mesh->GetSolid(), matPlasticScintillator, "PlasticScint_" + name + "_Logic");
    G4VPhysicalVolume* physPlasticScint2 = new G4PVPlacement(rotation, translation, logicPlasticScint2, "PlasticScint_" + name, logicWorld, false, 0, true);
    logicScintillators.push_back(logicPlasticScint2);
    
    // Cuvette
    auto cuvetteMesh1 = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_cuvetteParts_1_cuvette.stl");
    cuvetteMesh1->SetScale(1.0);
    G4LogicalVolume* logicCuvette1 = new G4LogicalVolume(cuvetteMesh1->GetSolid(), matCuvette, "Cuvette_" + name + "_Logic");
    G4VPhysicalVolume* physCuvette1 = new G4PVPlacement(rotation, translation, logicCuvette1, "Cuvette_" + name, logicWorld, false, 0, true);

    auto cuvetteMesh2 = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_cuvetteParts_2_cuvette.stl");
    cuvetteMesh2->SetScale(1.0);
    G4LogicalVolume* logicCuvette2 = new G4LogicalVolume(cuvetteMesh2->GetSolid(), matCuvette, "Cuvette_" + name + "_Logic");
    G4VPhysicalVolume* physCuvette2 = new G4PVPlacement(rotation, translation, logicCuvette2, "Cuvette_" + name, logicWorld, false, 0, true);

    // CuvettePVC
    auto cuvettePVCMesh1 = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_cuvetteParts_1_PVC.stl");
    cuvettePVCMesh1->SetScale(1.0);
    G4LogicalVolume* logicCuvettePVC1 = new G4LogicalVolume(cuvettePVCMesh1->GetSolid(), matPVC, "CuvettePVC_" + name + "_Logic");
    G4VPhysicalVolume* physCuvettePVC1 = new G4PVPlacement(rotation, translation, logicCuvettePVC1, "CuvettePVC_" + name, logicWorld, false, 0, true);

    auto cuvettePVCMesh2 = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_cuvetteParts_2_PVC.stl");
    cuvettePVCMesh2->SetScale(1.0);
    G4LogicalVolume* logicCuvettePVC2 = new G4LogicalVolume(cuvettePVCMesh2->GetSolid(), matPVC, "CuvettePVC_" + name + "_Logic");
    G4VPhysicalVolume* physCuvettePVC2 = new G4PVPlacement(rotation, translation, logicCuvettePVC2, "CuvettePVC_" + name, logicWorld, false, 0, true);

    // CuvetteLS
    auto cuvetteLSMesh1 = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_cuvetteParts_1_LS.stl");
    cuvetteLSMesh1->SetScale(1.0);
    G4LogicalVolume* logicCuvetteLS1 = new G4LogicalVolume(cuvetteLSMesh1->GetSolid(), matGenerator, "CuvetteLS_" + name + "_Logic");
    logicScintillators.push_back(logicCuvetteLS1);
    G4VPhysicalVolume* physCuvetteLS1 = new G4PVPlacement(rotation, translation, logicCuvetteLS1, "CuvetteLS_" + name, logicWorld, false, 0, true);

    auto cuvetteLSMesh2 = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_cuvetteParts_2_LS.stl");
    cuvetteLSMesh2->SetScale(1.0);
    G4LogicalVolume* logicCuvetteLS2 = new G4LogicalVolume(cuvetteLSMesh2->GetSolid(), matGenerator, "CuvetteLS_" + name + "_Logic");
    logicScintillators.push_back(logicCuvetteLS2);
    G4VPhysicalVolume* physCuvetteLS2 = new G4PVPlacement(rotation, translation, logicCuvetteLS2, "CuvetteLS_" + name, logicWorld, false, 0, true);

    // Ps_CoverDisk
    auto coverDiskMesh = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_ps_disk_cover_1_Body1.stl");
    coverDiskMesh->SetScale(1.0);
    G4LogicalVolume* logicCoverDisk = new G4LogicalVolume(coverDiskMesh->GetSolid(), matPhotopolymer, "CoverDisk_" + name + "_Logic");
    G4VPhysicalVolume* physCoverDisk = new G4PVPlacement(rotation, translation, logicCoverDisk, "CoverDisk_" + name, logicWorld, false, 0, true);

    auto holderDiskMesh = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_ps_disk_holder_1_Body1.stl");
    holderDiskMesh->SetScale(1.0);
    G4LogicalVolume* logicHolderDisk = new G4LogicalVolume(holderDiskMesh->GetSolid(), matPhotopolymer, "HolderDisk_" + name + "_Logic");
    G4VPhysicalVolume* physHolderDisk = new G4PVPlacement(rotation, translation, logicHolderDisk, "HolderDisk_" + name, logicWorld, false, 0, true);

    // SiPMs 1-4
    auto sipmMesh1 = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_SiPMs_1_Body1.stl");
    sipmMesh1->SetScale(1.0);
    G4LogicalVolume* logicSiPM_vol1 = new G4LogicalVolume(sipmMesh1->GetSolid(), matSiPM, "SiPMMount1_" + name + "_Logic");
    logicSiPM.push_back(logicSiPM_vol1);
    G4VPhysicalVolume* physSiPM1 = new G4PVPlacement(rotation, translation, logicSiPM_vol1, "SiPMMount1_" + name, logicWorld, false, 0, true);

    auto sipmMesh2 = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_SiPMs_1_Body2.stl");
    sipmMesh2->SetScale(1.0);
    G4LogicalVolume* logicSiPM_vol2 = new G4LogicalVolume(sipmMesh2->GetSolid(), matSiPM, "SiPMMount2_" + name + "_Logic");
    logicSiPM.push_back(logicSiPM_vol2);
    G4VPhysicalVolume* physSiPM2 = new G4PVPlacement(rotation, translation, logicSiPM_vol2, "SiPMMount2_" + name, logicWorld, false, 0, true);

    auto sipmMesh3 = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_SiPMs_1_Body3.stl");
    sipmMesh3->SetScale(1.0);
    G4LogicalVolume* logicSiPM_vol3 = new G4LogicalVolume(sipmMesh3->GetSolid(), matSiPM, "SiPMMount3_" + name + "_Logic");
    logicSiPM.push_back(logicSiPM_vol3);
    G4VPhysicalVolume* physSiPM3 = new G4PVPlacement(rotation, translation, logicSiPM_vol3, "SiPMMount3_" + name, logicWorld, false, 0, true);
    
    auto sipmMesh4 = CADMesh::TessellatedMesh::FromSTL(prefix + "Positronium_Generator_Positronium_Generator_SiPMs_1_Body5.stl");
    sipmMesh4->SetScale(1.0);
    G4LogicalVolume* logicSiPM_vol4 = new G4LogicalVolume(sipmMesh4->GetSolid(), matSiPM, "SiPMMount4_" + name + "_Logic");
    logicSiPM.push_back(logicSiPM_vol4);
    G4VPhysicalVolume* physSiPM4 = new G4PVPlacement(rotation, translation, logicSiPM_vol4, "SiPMMount4_" + name, logicWorld, false, 0, true);
    logicCalorimeter=logicSiPM_vol1; // Set the first SiPM as the reference for the sensitive detector
    // ==========================================
    // 7. OPTICAL SURFACES & BORDERS
    // ==========================================
    
    // We use 'static' to ensure these templates are only created ONCE in memory
    static G4OpticalSurface* surfTyvekWrap = nullptr;
    static G4OpticalSurface* surfSiPM = nullptr;
    static G4OpticalSurface* surfMetalFoil = nullptr;

    if (!surfTyvekWrap) {
        // --- 1 & 4. Reflector Setup (Tyvek on Ground Surface) ---
        // 'groundfrontpainted' models a rough surface (ground) wrapped in a reflector (Tyvek)
        surfTyvekWrap = new G4OpticalSurface("Tyvek_Surf");
        surfTyvekWrap->SetType(dielectric_metal); // Tyvek is opaque, treated as metal boundary for reflection
        surfTyvekWrap->SetModel(unified);
        surfTyvekWrap->SetFinish(groundfrontpainted); 
        surfTyvekWrap->SetSigmaAlpha(0.2); // Roughness parameter (can be tuned)

        G4MaterialPropertiesTable* mptTyvek = new G4MaterialPropertiesTable();
        // Setup energy array (adjust the range 2.0-3.5 eV based on your plastic scintillator emission spectrum)
        G4double photonEnergy[] = { 2.0*eV, 2.5*eV, 3.0*eV, 3.5*eV }; 
        G4double reflectivity[] = { 0.98, 0.98, 0.98, 0.98 }; // Tyvek typically has ~98% reflectivity
        mptTyvek->AddProperty("REFLECTIVITY", photonEnergy, reflectivity, 4);
        surfTyvekWrap->SetMaterialPropertiesTable(mptTyvek);
    }

    if (!surfSiPM) {
        // --- 3. SiPM Surface Setup ---
        // Ensure photons hitting the SiPM are absorbed (detected) and not reflected back
        surfSiPM = new G4OpticalSurface("SiPM_Surf");
        surfSiPM->SetType(dielectric_metal); 
        surfSiPM->SetModel(unified);
        surfSiPM->SetFinish(polished);
        
        G4MaterialPropertiesTable* mptSiPM = new G4MaterialPropertiesTable();
        G4double photonEnergy[] = { 2.0*eV, 2.5*eV, 3.0*eV, 3.5*eV }; 
        G4double reflectivity[] = { 0.0, 0.0, 0.0, 0.0 }; // 0% reflection = completely absorbed
        G4double efficiency[]   = { 1.0, 1.0, 1.0, 1.0 }; // 100% QE for now (can map to actual PDE later)
        mptSiPM->AddProperty("REFLECTIVITY", photonEnergy, reflectivity, 4);
        mptSiPM->AddProperty("EFFICIENCY", photonEnergy, efficiency, 4);
        surfSiPM->SetMaterialPropertiesTable(mptSiPM);
    }
    if (!surfMetalFoil) {
        surfMetalFoil = new G4OpticalSurface("MetalFoil_Surf");
        surfMetalFoil->SetType(dielectric_metal); // 金屬邊界
        surfMetalFoil->SetModel(unified);
        surfMetalFoil->SetFinish(polished);       // 假設為平滑金屬箔片
        
        G4MaterialPropertiesTable* mptMetalFoil = new G4MaterialPropertiesTable();
        G4double photonEnergy[] = { 2.0*eV, 2.5*eV, 3.0*eV, 3.5*eV }; 
        G4double reflectivity[] = { 0.50, 0.50, 0.50, 0.50 }; // 鈦金屬的可見光反射率約 50%
        mptMetalFoil->AddProperty("REFLECTIVITY", photonEnergy, reflectivity, 4);
        surfMetalFoil->SetMaterialPropertiesTable(mptMetalFoil);
    }


    // --- Apply Skin Surfaces ---
    // Apply Tyvek wrapping to the ground surfaces of the disks
    new G4LogicalSkinSurface("CoverDisk_Skin_" + name, logicCoverDisk, surfTyvekWrap);
    new G4LogicalSkinSurface("HolderDisk_Skin_" + name, logicHolderDisk, surfTyvekWrap);
    new G4LogicalSkinSurface("Cuvette1_Skin_" + name, logicCuvette1, surfTyvekWrap);
    new G4LogicalSkinSurface("Cuvette2_Skin_" + name, logicCuvette2, surfTyvekWrap);
    new G4LogicalSkinSurface("TiFoil_Skin_" + name, logicFoilsourceMesh, surfMetalFoil);

    // --- 2. Optical Borders (Air Gap Explanation) ---
    // NOTE: Since you explicitly stated there are air gaps between components (Air Coupled),
    // NO G4LogicalBorderSurface is required or allowed here. 
    // Geant4 automatically handles optical photon refraction & reflection at the boundaries 
    // between your physical volumes and the logicWorld (Air gap) using Snell's law.
    // 
    // [CRITICAL REMINDER]: Ensure your matAir (logicWorld material) and matPlasticScintillator 
    // both have the "RINDEX" (Refractive Index) array defined in their G4MaterialPropertiesTable 
    // in your material definition file, otherwise optical photons will instantly die at the boundaries!
}

void MyDetectorConstruction::ConstructCalorimeter() {
    // Place a single unit at origin
    if(is3DCalorimeter){
        // Generate custom coordinates
        std::ifstream coordFile(coordinate_name);
        if (!coordFile.is_open()) {
            G4cerr << "Error: Cannot open coordinates.txt for calorimeter positions!" << G4endl;
            return;  // Or fall back to old lattice code if preferred
        }

        std::string line;
        G4int counter = 0;
        while (std::getline(coordFile, line)) {
            std::istringstream linestream(line);
            // Skip empty lines or comments
            if (line.empty() || line[0] == '#') continue;

            std::istringstream iss(line);
            std::vector<G4double> vals;
            G4double val;
            while (iss >> val) {
                vals.push_back(val);
            }

            if (vals.size() != 3 && vals.size() != 6) {
                G4cout << "Warning: Skipping invalid line in coordinates.txt (expected 3 or 5 parameters): " << line << G4endl;
                continue;
            }

            G4double x = vals[0];
            G4double y = vals[1];
            G4double z = vals[2];
            G4double rot1 = 0.;
            G4double rot2 = 0.;
            G4double rot3= 0.;

            if (vals.size() == 6) {
                rot1 = vals[3];
                rot2 = vals[4];
                rot3=vals[5];
            }

            // Position in mm (adjust units if your file uses different, e.g., *cm)
            G4ThreeVector translation(x, y, z);

            // Unique name
            G4String name = "calor_unit_" + std::to_string(counter++);
            // Call the unit constructor with rotations
            ConstructCalorimeter_unit_3d(translation,  name,rot1 * deg, rot2 * deg, rot3 * deg);
        }

        coordFile.close();
    }
    else{
        ConstructCalorimeter_unit(G4ThreeVector(0,0,0), 0*deg, "");
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

void MyDetectorConstruction::SetLightYield(G4double val) {
    fLightYield = val;
    G4cout << "Updating Liquid Scintillator Light Yield to: " << fLightYield << " photons/MeV" << G4endl;
    
    // If the material table exists, update the specific property
    if (matWater && matWater->GetMaterialPropertiesTable()) {
        matWater->GetMaterialPropertiesTable()->AddConstProperty("SCINTILLATIONYIELD", fLightYield / MeV);
        
        // Inform the RunManager that cross-sections/physics must be recalculated
        G4RunManager::GetRunManager()->PhysicsHasBeenModified(); 
    }
}