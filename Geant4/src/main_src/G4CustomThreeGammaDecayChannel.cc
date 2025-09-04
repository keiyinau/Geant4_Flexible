#include "G4CustomThreeGammaDecayChannel.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4LorentzVector.hh"
#include <cmath>
#include <random>

G4CustomThreeGammaDecayChannel::G4CustomThreeGammaDecayChannel(const G4String& theParentName, G4double theBR)
    : G4VDecayChannel("Custom 3 Gamma Decay", theParentName, theBR, 3, "gamma", "gamma", "gamma") {}

G4DecayProducts* G4CustomThreeGammaDecayChannel::DecayIt(G4double parentMass) {
    if (parentMass <= 0.0) parentMass = GetParentMass();  // Fallback to default if not provided

    G4double M = parentMass;  // o-Ps rest mass (~1.022 MeV)
    auto dXS=[M](double e1,double e2, double e3){const double pi=std::acos(-1.0);
    return (((M/2 - e1) / (e2 * e3)) * ((M/2 - e1) / (e2 * e3)) +
                ((M/2 - e2) / (e1 * e3)) * ((M/2 - e2) / (e1 * e3)) +
                ((M/2 - e3) / (e1 * e2)) * ((M/2 - e3) / (e1 * e2))) /
               (pi * pi - 9);
            };


    // Sample energies with rejection for phase space
    G4double E1, E2, E3, p1, p2, p3, s;
    G4double rd1, rd2,rd3;
    G4float dXS_val, dXS_min, dXS_max, rng_gen;
    dXS_min=6.620097e-6;
    dXS_max=8.807577e-06;
    s=0;
    std::random_device rd;
    std::mt19937 gen(rd());
    while(s==0){
        rd1 = G4UniformRand();
        rd2 = G4UniformRand();
        rd3 = G4UniformRand();
        E1=0+rd1*M/2;
        E2=M/2-E1+rd2*E1;
        E3=M-E1-E2;
        if(E1>0 & E2>0 & E3>0){
            if(E1<M/2 & E2<M/2 & E3<M/2 & E1+E2+E3==M){
                rng_gen=dXS_min+(dXS_max-dXS_min)*rd3;
                dXS_val=dXS(E1,E2,E3);
                if(rng_gen<dXS_val){
                    std::vector<double> E_pshuffle = {E1, E2, E3};
                    std::shuffle(E_pshuffle.begin(), E_pshuffle.end(), gen);
                    E1=E_pshuffle[0];
                    E2=E_pshuffle[1];
                    E3=E_pshuffle[2];
                    s+=1;
                }
            }
        }
    }
        // Create parent at rest
    G4DynamicParticle* parentParticle = new G4DynamicParticle(GetParent(), G4ThreeVector(0., 0., 0.), 0.0);
    G4DecayProducts* products = new G4DecayProducts(*parentParticle);
    delete parentParticle;
    // Sample directions (photon 1, 2, 3 correspond to E1, E2, E3)
    // Random isotropic direction for photon 1
    G4double cost = 2.0 * G4UniformRand() - 1.0;
    G4double sint = std::sqrt(1.0 - cost * cost);
    G4double phi = twopi * G4UniformRand();
    G4ThreeVector dir1(sint * std::cos(phi), sint * std::sin(phi), cost);
    G4ThreeVector p1v = E1 * dir1;

    // Fixed cos(theta12) for angle between photon 1 and 2
    G4double cos_theta12 = (E3 * E3 - E1 * E1 - E2 * E2) / (2.0 * E1 * E2);
    // Ensure cos_theta12 is in valid range [-1, 1]
    if (std::abs(cos_theta12) > 1.0) {
        cos_theta12 = (cos_theta12 > 1.0) ? 1.0 : -1.0;
    }
    G4double theta12 = std::acos(cos_theta12);
    G4double sin_theta12 = std::sin(theta12);

    // Random azimuthal angle for photon 2 relative to photon 1
    G4double phi2 = twopi * G4UniformRand();

    // Construct dir2 in frame where dir1 is along z
    G4ThreeVector dir2(sin_theta12 * std::cos(phi2), sin_theta12 * std::sin(phi2), cos_theta12);

    // Rotate dir2 to align with actual dir1 (new z-axis)
    dir2.rotateUz(dir1);
    G4ThreeVector p2v = E2 * dir2;

    // Photon 3 momentum (conservation)
    G4ThreeVector p3v = -p1v - p2v;

    // Create photons
    G4DynamicParticle* gamma1 = new G4DynamicParticle(G4Gamma::Gamma(), p1v);
    G4DynamicParticle* gamma2 = new G4DynamicParticle(G4Gamma::Gamma(), p2v);
    G4DynamicParticle* gamma3 = new G4DynamicParticle(G4Gamma::Gamma(), p3v);

    products->PushProducts(gamma1);
    products->PushProducts(gamma2);
    products->PushProducts(gamma3);

#ifdef G4VERBOSE
    if (GetVerboseLevel() > 1) {
        G4cout << "G4CustomThreeGammaDecayChannel::DecayIt - Created 3 photons with energies: "
               << E1 / keV << ", " << E2 / keV << ", " << E3 / keV << " keV" << G4endl;
    }
#endif

    return products;
}