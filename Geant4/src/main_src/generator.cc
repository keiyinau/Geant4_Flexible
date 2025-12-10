#include "generator.hh"

MyPrimaryGenerator::MyPrimaryGenerator(){
	ParticleDefinition();

	// G4ParticleGun(number of particles for 1 event), here set 1 particle for 1 event
	fParticleGun = new G4ParticleGun(1);
	fParticleSource = new G4GeneralParticleSource();

	// Define fIon to be Ground state Na-22 (Z=11, A=22)
	// Co-60 (Z=27, A=60) can also be used for testing
	// Cs-137 (Z=55, A=137) can also be used for testing
	
	Z = 11;					// Atomic number (Proton Number)
	A = 22;					// Mass number
	ex_energy = 0.*keV;		// Excitation energy

	// Select particle generator. Options: 0 = fParticleSource, 1 = fParticleGun
	PS_or_PG = 0;

	// Set the default parameters for the fParticleSource
	pdParticleSource = fGeantino;								//options: fGamma, fPositron, fGeantino = fIon
	chargeParticleSource = 0.*eplus;
	fParticleSource->SetParticleCharge(chargeParticleSource);
	fParticleSource->SetParticleDefinition(pdParticleSource);

	// Set the default parameters for the fParticleGun
	pdParticleGun = fGamma;									//options: fGamma, fPositron, fGeantino = fIon, fo_Ps, fp_Ps
	//posParticleGun = G4ThreeVector(0.*cm, 0.*cm, 0.*cm);
	posParticleGun = G4ThreeVector(0.*cm, 0.*cm, 0.*cm);
	momDirectionParticleGun = G4ThreeVector(0., 0., 1.);
	kinParticleGun = 100*keV; 
	chargeParticleGun = 0.*eplus;
	//if (pdParticleGun == fPositron){
	//	G4double maxEnergy = 600. * keV;  // Na-22 beta+ endpoint
	//	G4double energy = maxEnergy * G4UniformRand();  // Simple uniform for demo; use Fermi function for accurate spectrum
//
	//	// Direction: isotropic or along a beam
	//	G4double theta = std::acos(2 * G4UniformRand() - 1);
	//	G4double phi = 2 * CLHEP::pi * G4UniformRand();
	//	G4ThreeVector direction(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
//
	//	fParticleGun->SetParticleDefinition(pdParticleGun);
	//	fParticleGun->SetParticleEnergy(energy);
	//	fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));  // At source position
	//	fParticleGun->SetParticleMomentumDirection(direction);
//
	//	// Set longitudinal polarization: +1 along momentum (for positrons)
	//	G4ThreeVector pol = direction.unit();  // Helicity +1
	//	fParticleGun->SetParticlePolarization(pol);
//
	//	fParticleGun->GeneratePrimaryVertex(anEvent);
	//}
	fParticleGun->SetParticlePosition(posParticleGun);
	fParticleGun->SetParticleMomentumDirection(momDirectionParticleGun);
	//fParticleGun->SetParticleMomentum(momParticleGun);
	fParticleGun->SetParticleEnergy(kinParticleGun);			//set the KineticEnergy of particle
	fParticleGun->SetParticleCharge(chargeParticleGun);
	fParticleGun->SetParticleDefinition(pdParticleGun);
}

MyPrimaryGenerator::~MyPrimaryGenerator(){
	delete fParticleGun;
	delete fParticleSource;
}

// Define Particle Definition
void MyPrimaryGenerator::ParticleDefinition(){
	fGamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
	fPositron = G4ParticleTable::GetParticleTable()->FindParticle("e+");
	fElectron=G4ParticleTable::GetParticleTable()->FindParticle("e-");
	fGeantino = G4ParticleTable::GetParticleTable()->FindParticle("geantino");
	fOpticalPhoton = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
	fo_Ps = G4OrthoPositronium::Definition();
	fp_Ps = G4ParaPositronium::Definition();
}

//This function is called at the begining of every event
void MyPrimaryGenerator::GeneratePrimaries(G4Event* anEvent){
	G4ParticleDefinition* pd_ParticleSource = fParticleSource->GetParticleDefinition();
	G4ParticleDefinition* pd_ParticleGun = fParticleGun->GetParticleDefinition();

	if (PS_or_PG == 0) {
		if (pd_ParticleSource == G4Geantino::Geantino()) {
			fIon = G4IonTable::GetIonTable()->GetIon(Z, A, ex_energy);
			fParticleSource->SetParticleDefinition(fIon);

			//fParticleSource->SetVerbosity(2);				// Turn on verbose mode and set the level 2, for debugging
			fParticleSource->GeneratePrimaryVertex(anEvent);
		}else
			fParticleSource->GeneratePrimaryVertex(anEvent);
	}
	else {
		if (pd_ParticleGun == G4Geantino::Geantino()) {
			fIon = G4IonTable::GetIonTable()->GetIon(Z, A, ex_energy);
			fParticleGun->SetParticleDefinition(fIon);
			fParticleGun->GeneratePrimaryVertex(anEvent);
		}else
			fParticleGun->GeneratePrimaryVertex(anEvent);
	}
}
