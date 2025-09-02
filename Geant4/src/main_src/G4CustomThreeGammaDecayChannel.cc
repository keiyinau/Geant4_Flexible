#include "G4CustomThreeGammaDecayChannel.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4LorentzVector.hh"
#include <cmath>

G4CustomThreeGammaDecayChannel::G4CustomThreeGammaDecayChannel(const G4String& theParentName, G4double theBR)
    : G4VDecayChannel("Custom 3 Gamma Decay", theParentName, theBR, 3, "gamma", "gamma", "gamma") {}

G4DecayProducts* G4CustomThreeGammaDecayChannel::DecayIt(G4double parentMass) {
    if (parentMass <= 0.0) parentMass = GetParentMass();  // Fallback to default if not provided

    G4double M = parentMass;  // o-Ps rest mass (~1.022 MeV)

    // Sample energies with rejection for phase space
    G4double E1, E2, E3, p1, p2, p3, maxp, sump;
    G4double rd1, rd2;
    do {
        rd1 = G4UniformRand();
        rd2 = G4UniformRand();
        if (rd2 > rd1) { G4double temp = rd1; rd1 = rd2; rd2 = temp; }  // Ensure rd2 <= rd1
        E1 = rd2 * M;
        E2 = (rd1 - rd2) * M;
        E3 = (1.0 - rd1) * M;
        p1 = E1; p2 = E2; p3 = E3;  // Massless: p = E
        maxp = std::max(E1, std::max(E2, E3));
        sump = E1 + E2 + E3;  // Should be ~M
    } while (maxp > sump - maxp);  // Equivalent to max(Ei) > M/2

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