#include "construction.hh"
#include "CADMesh.hh"
MyDetectorConstruction::MyDetectorConstruction() {
	// Define required materials
    logicOptical=true;
	DefineMaterials();


    coordinate_name="coordinates.txt";
    fLightYield=10000;		// Light yield of scintillator, in photons/MeV. This is a made-up value for demonstration; adjust based on actual material properties.
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
	//Define the world material as Air
	Air = nist->FindOrBuildMaterial("G4_AIR");
    std::vector<G4double> Air_absorption_Energy, Air_absorption_Index;
    readAndProcessData_Energy_cm_txt("AbsorptionLength_Air.txt", Air_absorption_Energy, Air_absorption_Index);
    G4MaterialPropertiesTable* mptAir = new G4MaterialPropertiesTable();
    mptAir->AddProperty("RINDEX", "Air");
    mptAir->AddProperty("ABSLENGTH", Air_absorption_Energy, Air_absorption_Index,Air_absorption_Index.size());
    

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
    // === NEW: Add Liquid Scintillation Properties ===
    // 1. Placeholder Emission Spectrum (e.g., 2.5 eV to 3.0 eV -> ~413nm to ~496nm)
    G4double liquidEnergy[] = { 2.5*eV, 3.0*eV }; 
    G4double liquidEmission[] = { 1.0, 1.0 }; // Flat emission profile
    mptWater->AddProperty("SCINTILLATIONCOMPONENT1", liquidEnergy, liquidEmission, 2);

    // 2. The Variable you will study: Light Yield (Photons per MeV)
    mptWater->AddConstProperty("SCINTILLATIONYIELD", fLightYield/MeV);
    
    // 3. Resolution and Timing (Typical liquid scintillator timings)
    mptWater->AddConstProperty("RESOLUTIONSCALE", 1.0);
    mptWater->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 5.0 * ns); // Fast decay
    
    matWater->SetMaterialPropertiesTable(mptWater);
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
    //mptLSO->AddConstProperty("RESOLUTIONSCALE", 0.);
    //mptLSO->AddProperty("SCINTILLATIONCOMPONENT1", LSO_emission_Energy, LSO_emission_fractions,LSO_emission_fractions.size());
    mptLSO->AddProperty("RINDEX", LSO_refraction_Energy, LSO_refraction_Index,LSO_refraction_Index.size());
    mptLSO->AddProperty("ABSLENGTH", LSO_absorption_Energy, LSO_absorption_Index,LSO_absorption_Index.size());
    //mptLSO->AddConstProperty("SCINTILLATIONYIELD", 26/keV);
    //mptLSO->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 40.0*ns);
    // End LSO

    // Define elements (use NIST for common ones)
    G4Element* elLu = nist->FindOrBuildElement("Lu");
    G4Element* elY  = nist->FindOrBuildElement("Y");
    G4Element* elSi = nist->FindOrBuildElement("Si");
    G4Element* elO  = nist->FindOrBuildElement("O");

    // Optionally include Ce for doping (trace amount; adjust fraction if needed)
    G4Element* elCe =  nist->FindOrBuildElement("Ce");

    // Define the material with density ~7.1 g/cm³ (common value for LYSO(Ce))
    matLYSO = new G4Material("LYSO", 7.1 * g/cm3, 4);  // 4 components (or 5 if including Ce)

    // Add elements by mass fraction (calculated from atomic ratios)
    matLYSO->AddElement(elLu, 0.7146);  // ~71.46%
    matLYSO->AddElement(elY, 0.0403);   // ~4.03%
    matLYSO->AddElement(elSi, 0.0637);  // ~6.37%
    matLYSO->AddElement(elO, 0.1814);   // ~18.14%
    G4MaterialPropertiesTable* mptLYSO = new G4MaterialPropertiesTable();
    std::vector<G4double> LYSO_emission_Energy, LYSO_emission_fractions;
    readAndProcessData_Energy("EmissionSpectrum_LYSO_Ce.csv", LYSO_emission_Energy, LYSO_emission_fractions);
    std::vector<G4double> LYSO_refraction_Energy, LYSO_refraction_Index;
    readAndProcessData_txt("RefractiveIndex_LYSO_Ce.txt", LYSO_refraction_Energy, LYSO_refraction_Index);
    std::vector<G4double> LYSO_absorption_Energy, LYSO_absorption_Index;
    readAndProcessData_Energy_cm_txt("AbsorptionLength_LYSO_Ce.txt", LYSO_absorption_Energy, LYSO_absorption_Index);
    std::vector<G4double> LYSO_LY_Nonproportion_Energy, LYSO_LY_Nonproportion_relative;
    readAndProcessData_Nonproportionality("Nonproportionality_LYSO_Ce_Relative.txt", LYSO_LY_Nonproportion_Energy, LYSO_LY_Nonproportion_relative);
    G4double baseYield=33.0000/keV; //Previous 33 keV
    std::vector<G4double> LYSO_LY_Nonproportion_fractions(LYSO_LY_Nonproportion_relative.size());
    for(int i=0;i<LYSO_LY_Nonproportion_relative.size();i++){
        LYSO_LY_Nonproportion_fractions[i]=LYSO_LY_Nonproportion_relative[i]*baseYield;
    }

    mptLYSO->AddConstProperty("SCINTILLATIONYIELD", baseYield); 
    
    mptLYSO->AddProperty("ELECTRONSCINTILLATIONYIELD", LYSO_LY_Nonproportion_Energy, LYSO_LY_Nonproportion_fractions, LYSO_LY_Nonproportion_fractions.size());
    mptLYSO->AddConstProperty("ELECTRONSCINTILLATIONYIELD1", 1.0);
    mptLYSO->AddConstProperty("RESOLUTIONSCALE", 0);
    mptLYSO->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 40. * ns);
    mptLYSO->AddProperty("SCINTILLATIONCOMPONENT1", LYSO_emission_Energy, LYSO_emission_fractions,LYSO_emission_fractions.size());
    mptLYSO->AddProperty("RINDEX", LYSO_refraction_Energy, LYSO_refraction_Index,LYSO_refraction_Index.size());
    mptLYSO->AddProperty("ABSLENGTH", LYSO_absorption_Energy, LYSO_absorption_Index,LYSO_absorption_Index.size());
    //mptLYSO->AddConstProperty("BIRKS_ETA_H", 0.002,true);
    //mptLYSO->AddConstProperty("ONSAGER_ETA_EH", 0.81,true);
    //mptLYSO->AddConstProperty("ONSAGER_k_O", 0.29,true);
    //matLYSO->GetIonisation()->SetBirksConstant(0.186 * mm/MeV); // For birks-onsager
    //matLYSO->GetIonisation()->SetBirksConstant(0.0028 * cm/MeV);

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
	// End CsI

	// Define Aluminium for wrapping and protection
	matAl = nist->FindOrBuildMaterial("G4_Al");
	G4MaterialPropertiesTable* mptAl = new G4MaterialPropertiesTable();
    const G4int nEntries = 2; 
    G4double PhotonEnergy[nEntries] = {1.5 * eV, 3.0 * eV}; 
    G4double RIndex_al[nEntries] = {1.37, 0.44}; 
    mptAl->AddProperty("RINDEX", PhotonEnergy, RIndex_al, nEntries);
    // End Aluminium

	// Define Acrylic
	matAcrylic = nist->FindOrBuildMaterial("G4_PLEXIGLASS");
	G4MaterialPropertiesTable* mptAcrylic = new G4MaterialPropertiesTable();
    mptAcrylic->AddProperty("RINDEX", "PMMA");
	// End Acrylic
    
    // Define Glass
    matGlass = nist->FindOrBuildMaterial("G4_Pyrex_Glass");
    //// 2. Define photon energies (e.g., covering the visible spectrum)
    G4double photonEnergies[] = { 2.034*eV, 2.384*eV, 2.755*eV, 3.100*eV }; // ~610nm to 400nm
    const G4int numEntries = sizeof(photonEnergies)/sizeof(G4double);

    // 3. Define the Refractive Index of Pyrex/Borosilicate glass (~1.47)
    G4double rIndexGlass[] = { 1.47, 1.47, 1.47, 1.47 }; 

    // 4. Create the Material Properties Table and attach it
    G4MaterialPropertiesTable* mptGlass = new G4MaterialPropertiesTable();
    mptGlass->AddProperty("RINDEX", photonEnergies, rIndexGlass, numEntries);
    // 5. Attach the table to the NIST material
    matGlass->SetMaterialPropertiesTable(mptGlass);
    // End Glass

	// Define Tapflon(teflon) for wrapping
	matTeflon = nist->FindOrBuildMaterial("G4_TEFLON");
	std::vector<G4double> tapflon_reflectance_Energy, tapflon_reflectance_fractions;
	readAndProcessData_Energy_txt("teflon_Reflectance-modified.txt", tapflon_reflectance_Energy, tapflon_reflectance_fractions);
    std::vector<G4double> tapflon_refraction_Energy, tapflon_refraction_Index;
	readAndProcessData_Energy_txt("Refraction_Index_Teflon_Gray.txt", tapflon_refraction_Energy, tapflon_refraction_Index);
	mptTeflon = new G4MaterialPropertiesTable();
	mptTeflon->AddProperty("REFLECTIVITY", tapflon_reflectance_Energy, tapflon_reflectance_fractions,tapflon_reflectance_fractions.size());
    mptTeflon->AddProperty("RINDEX", tapflon_refraction_Energy, tapflon_refraction_Index,tapflon_refraction_Index.size());
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
	//mptSi->AddProperty("REFLECTIVITY", Si_reflectance_Energy, Si_reflectance_fractions,Si_reflectance_fractions.size());
    mptSi->AddProperty("TRANSMITTANCE", Si_transmission_Energy, Si_rtransmission_Index,Si_rtransmission_Index.size());	
    mptSi->AddProperty("RINDEX", Si_refraction_Energy, Si_refraction_Index,Si_refraction_Energy.size());	

    // CsI-Teflon (reflective surface)
    //surfCsI_Teflon = new G4OpticalSurface("CsI_Teflon_Surface");
    ////surfCsI_Teflon->SetType(dielectric_dielectric); // Teflon as reflective surface, dielectric_metal
    ////surfCsI_Teflon->SetModel(unified);
    ////surfCsI_Teflon->SetFinish(polished);
    //surfCsI_Teflon->SetType(dielectric_dielectric);
    //surfCsI_Teflon->SetModel(unified);
    //surfCsI_Teflon->SetFinish(ground); // specular->polishedteflonair, diffusive->groundteflonair
    //surfCsI_Teflon->SetSigmaAlpha(0.2);
    ////End surface
//
    //// CsI-SiPM (dielectric-dielectric interface)
    //surfCsI_SiPM = new G4OpticalSurface("CsI_SiPM_Surface");
    //surfCsI_SiPM->SetType(dielectric_dielectric);
    //surfCsI_SiPM->SetModel(glisur); // Glisur for smooth dielectric interface
    //surfCsI_SiPM->SetFinish(polished);
	//// End SiPM
//
    //// CsI-AlFoil (reflective surface)
    //surfCsI_AlFoil = new G4OpticalSurface("CsI_AlFoil_Surface");
    //surfCsI_AlFoil->SetType(dielectric_metal); // Al as reflective surface
    //surfCsI_AlFoil->SetModel(unified);
    //surfCsI_AlFoil->SetFinish(ground);


    G4OpticalSurface* surfCrystalGrease = new G4OpticalSurface("CrystalGrease");
    surfCrystalGrease->SetType(dielectric_dielectric);
    surfCrystalGrease->SetFinish(polished);
    surfCrystalGrease->SetModel(unified);
    surfCrystalGrease->SetSigmaAlpha(0.05*degree);
    //End surface
    if(logicOptical){
        Air->SetMaterialPropertiesTable(mptAir);
        matWater->SetMaterialPropertiesTable(mptWater);
        matLSO->SetMaterialPropertiesTable(mptLSO);
        matLYSO->SetMaterialPropertiesTable(mptLYSO);
        matSi->SetMaterialPropertiesTable(mptSi);
        matCsI->SetMaterialPropertiesTable(mptCsI);
        matAl->SetMaterialPropertiesTable(mptAl);
        matAcrylic->SetMaterialPropertiesTable(mptAcrylic);
        matTeflon->SetMaterialPropertiesTable(mptTeflon);
    }


    std::cout<<"==========================="<<std::endl;
    std::cout<<"Printing the material properties of LYSO"<<std::endl;
    mptLYSO->DumpTable();
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
    
    std::string prefix = "pixelated_lyso/";
    
    // === 1. MAIN LYSO CRYSTAL ===
    auto mainLysoMesh = CADMesh::TessellatedMesh::FromSTL(prefix + "Pixelized_Test_Pixelized_Test_LYSO_twoend_1_LYSO.stl");
    mainLysoMesh->SetScale(1.0);
    G4LogicalVolume* logicMainLYSO = new G4LogicalVolume(mainLysoMesh->GetSolid(), matLYSO, "MainLYSO_" + name + "_Logic");
    logicScintillators.push_back(logicMainLYSO); 
    logicCalorimeter=logicMainLYSO;
    // [FIXED] Added G4VPhysicalVolume* to store the pointer
    G4VPhysicalVolume* physMainLYSO = new G4PVPlacement(rotation, translation, logicMainLYSO, "MainLYSO_" + name, logicWorld, false, 0, true);
    
    auto mainWrapMesh = CADMesh::TessellatedMesh::FromSTL(prefix + "Pixelized_Test_Pixelized_Test_LYSO_twoend_1_Alwrap.stl");
    mainWrapMesh->SetScale(1.0);
    G4LogicalVolume* logicMainWrap = new G4LogicalVolume(mainWrapMesh->GetSolid(), matWrapping, "MainWrap_" + name + "_Logic");
    logicTapflon.push_back(logicMainWrap); 
    // [FIXED] Added G4VPhysicalVolume* to store the pointer
    G4VPhysicalVolume* physMainWrap = new G4PVPlacement(rotation, translation, logicMainWrap, "MainWrap_" + name, logicWorld, false, 0, true);

    auto mainSiPM1 = CADMesh::TessellatedMesh::FromSTL(prefix + "Pixelized_Test_Pixelized_Test_LYSO_twoend_1_SiPM1.stl");
    mainSiPM1->SetScale(1.0);
    G4LogicalVolume* logicMainSiPM1 = new G4LogicalVolume(mainSiPM1->GetSolid(), matSiPM, "MainSiPM1_" + name + "_Logic");
    logicSiPM.push_back(logicMainSiPM1); 
    // [FIXED] Added G4VPhysicalVolume* to store the pointer
    G4VPhysicalVolume* physMainSiPM1 = new G4PVPlacement(rotation, translation, logicMainSiPM1, "MainSiPM1_" + name, logicWorld, false, 0, true);

    auto mainSiPM2 = CADMesh::TessellatedMesh::FromSTL(prefix + "Pixelized_Test_Pixelized_Test_LYSO_twoend_1_SiPM2.stl");
    mainSiPM2->SetScale(1.0);
    G4LogicalVolume* logicMainSiPM2 = new G4LogicalVolume(mainSiPM2->GetSolid(), matSiPM, "MainSiPM2_" + name + "_Logic");
    logicSiPM.push_back(logicMainSiPM2); 
    // [FIXED] Added G4VPhysicalVolume* to store the pointer
    G4VPhysicalVolume* physMainSiPM2 = new G4PVPlacement(rotation, translation, logicMainSiPM2, "MainSiPM2_" + name, logicWorld, false, 0, true);

    // === 3. ONLY THE MIDDLE SIDE LYSO CRYSTAL ===
    auto sideLysoMesh = CADMesh::TessellatedMesh::FromSTL(prefix + "Pixelized_Test_Pixelized_Test_LYSO_single_1_LYSO.stl");
    sideLysoMesh->SetScale(1.0);
    G4LogicalVolume* logicSideLYSO = new G4LogicalVolume(sideLysoMesh->GetSolid(), matLSO, "SideLYSO" + name + "_Logic");
    // [FIXED] Added G4VPhysicalVolume* to store the pointer
    G4VPhysicalVolume* physSideLYSO = new G4PVPlacement(rotation, translation, logicSideLYSO, "SideLYSO" + name, logicWorld, false, 0, true);

    auto wrapMesh = CADMesh::TessellatedMesh::FromSTL(prefix + "Pixelized_Test_Pixelized_Test_LYSO_single_1_wrapping.stl");
    wrapMesh->SetScale(1.0);
    G4LogicalVolume* logicWrap = new G4LogicalVolume(wrapMesh->GetSolid(), matWrapping, "WrapSide" + name + "_Logic");
    logicTapflon.push_back(logicWrap); 
    // [FIXED] Added G4VPhysicalVolume* to store the pointer
    G4VPhysicalVolume* physSideWrap = new G4PVPlacement(rotation, translation, logicWrap, "WrapSide" + name, logicWorld, false, 0, true);

    auto sideSiPM1 = CADMesh::TessellatedMesh::FromSTL(prefix + "Pixelized_Test_Pixelized_Test_LYSO_single_1_SiPM.stl");
    sideSiPM1->SetScale(1.0);
    G4LogicalVolume* logicSideSiPM1 = new G4LogicalVolume(sideSiPM1->GetSolid(), matSiPM, "SideSiPM1" + name + "_Logic");
    logicSiPM.push_back(logicSideSiPM1); 
    // [FIXED] Added G4VPhysicalVolume* to store the pointer
    G4VPhysicalVolume* physSideSiPM1 = new G4PVPlacement(rotation, translation, logicSideSiPM1, "SideSiPM1" + name, logicWorld, false, 0, true);

    // ==========================================
    // OPTICAL SURFACES
    // ==========================================
    
    // 1. 定義 LYSO 與 Wrapping (Al Foil) 之間的表面
    G4OpticalSurface* wrapOpSurface = new G4OpticalSurface("WrapSurface_" + name);
    wrapOpSurface->SetType(dielectric_metal); 
    wrapOpSurface->SetModel(unified);
    wrapOpSurface->SetFinish(ground); 
    wrapOpSurface->SetSigmaAlpha(0.1); 

    G4MaterialPropertiesTable* mptWrapSurface = new G4MaterialPropertiesTable();
    const G4int numEntries = 2;
    G4double photonEnergy[numEntries] = {1.5 * eV, 3.0 * eV}; 
    G4double reflectivity[numEntries] = {0.90, 0.90}; 

    mptWrapSurface->AddProperty("REFLECTIVITY", photonEnergy, reflectivity, numEntries);
    wrapOpSurface->SetMaterialPropertiesTable(mptWrapSurface);

    // 2. 定義 LYSO 與 SiPM 之間的表面 (假設中間有光學硅脂或完美接觸)
    // [FIXED] 先 new 出 sipmOpSurface，然後再設定 Grease 屬性
    G4OpticalSurface* sipmOpSurface = new G4OpticalSurface("SiPMSurface_" + name);
    sipmOpSurface->SetType(dielectric_dielectric);
    sipmOpSurface->SetModel(unified);
    sipmOpSurface->SetFinish(polished); 
    

    G4MaterialPropertiesTable* mptGrease = new G4MaterialPropertiesTable();
    const G4int nEntriesGrease = 2;
    G4double energyGrease[nEntriesGrease] = {1.5 * eV, 3.0 * eV}; 
    G4double rindexGrease[nEntriesGrease] = {1.46, 1.46}; 

    mptGrease->AddProperty("RINDEX", energyGrease, rindexGrease, nEntriesGrease);
    
    // 將 Grease 的屬性附加到剛才建立的 Surface 上
    sipmOpSurface->SetMaterialPropertiesTable(mptGrease);

    // 3. 建立 Logical Border Surfaces
    new G4LogicalBorderSurface("MainLYSO_Wrap_Surf_" + name, physMainLYSO, physMainWrap, wrapOpSurface);
    new G4LogicalBorderSurface("SideLYSO_Wrap_Surf_" + name, physSideLYSO, physSideWrap, wrapOpSurface);

    new G4LogicalBorderSurface("MainLYSO_SiPM1_Surf_" + name, physMainLYSO, physMainSiPM1, sipmOpSurface);
    new G4LogicalBorderSurface("MainLYSO_SiPM2_Surf_" + name, physMainLYSO, physMainSiPM2, sipmOpSurface);
    new G4LogicalBorderSurface("SideLYSO_SiPM1_Surf_" + name, physSideLYSO, physSideSiPM1, sipmOpSurface);
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