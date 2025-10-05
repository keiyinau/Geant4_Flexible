#include "G4CustomThreeGammaDecayChannel.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4LorentzVector.hh"
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>

G4CustomThreeGammaDecayChannel::G4CustomThreeGammaDecayChannel(const G4String& theParentName, G4double theBR)
    : G4VDecayChannel("Custom 3 Gamma Decay", theParentName, theBR, 3, "gamma", "gamma", "gamma") {}

G4DecayProducts* G4CustomThreeGammaDecayChannel::DecayIt(G4double parentMass) {
    if (parentMass <= 0.0) parentMass = GetParentMass();  // Fallback to default if not provided

    G4double M = parentMass;  // o-Ps rest mass (~1.022 MeV)
    const static G4double min_energy = 0.001 * M / 2;  // Infrared cutoff, ~0.5 keV normalized to M/2
    const G4int n_points = 1;  // Generate one event per call
    const G4int n_probe = 1000000;  // Number of probes to estimate max density

    // Dalitz density function
    auto dalitz_density = [M](G4double e1, G4double e2, G4double e3) {
        G4double z = M - e1 - e2;
        if (e1 < min_energy || e2 < min_energy || z < min_energy) return 0.0;
        if (e1 >= M / 2 || e2 >= M / 2 || z >= M / 2) return 0.0;
        G4double term1 = std::pow((M / 2 - e1), 2) / (std::pow(e2, 2) * std::pow(z, 2));
        G4double term2 = std::pow((M / 2 - e2), 2) / (std::pow(e1, 2) * std::pow(z, 2));
        G4double term3 = std::pow((M / 2 - z), 2) / (std::pow(e1, 2) * std::pow(e2, 2));
        return term1 + term2 + term3;
    };

    // Find approximate max density for rejection sampling
    G4double max_d = 0.0;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (G4int i = 0; i < n_probe; ++i) {
        G4double x = min_energy + dis(gen) * (M / 2 - 2 * min_energy);
        G4double y = min_energy + dis(gen) * (M / 2 - 2 * min_energy);
        G4double z = M - x - y;
        if (z > min_energy && z < M / 2 && x + y < M) {
            G4double d = dalitz_density(x, y, z);
            if (d > max_d) max_d = d;
        }
    }
    max_d *= 1.1;  // Safety margin

    // Generate one event using rejection sampling
    G4double E1, E2, E3;
    G4int accepted = 0;
    while (accepted == 0) {
        G4double x = min_energy + dis(gen) * (M / 2 - min_energy);
        G4double y = min_energy + dis(gen) * (M / 2 - min_energy);
        E3 = M - x - y;
        if (E3 > min_energy && E3 < M / 2 && x + y < M) {
            G4double d = dalitz_density(x, y, E3);
            if (dis(gen) * max_d < d) {
                // Shuffle energies to avoid bias
                std::vector<G4double> E_pshuffle = {x, y, E3};
                std::shuffle(E_pshuffle.begin(), E_pshuffle.end(), gen);
                E1 = E_pshuffle[0];
                E2 = E_pshuffle[1];
                E3 = E_pshuffle[2];
                accepted = 1;
            }
        }
    }

    // Create parent at rest
    G4DynamicParticle* parentParticle = new G4DynamicParticle(GetParent(), G4ThreeVector(0., 0., 0.), 0.0);
    G4DecayProducts* products = new G4DecayProducts(*parentParticle);
    delete parentParticle;

    // Sample directions
    // Photon 1 direction (isotropic)
    G4double cost = 2.0 * G4UniformRand() - 1.0;
    G4double sint = std::sqrt(1.0 - cost * cost);
    G4double phi = 2.0 * M_PI * G4UniformRand();
    G4ThreeVector dir1(sint * std::cos(phi), sint * std::sin(phi), cost);
    G4ThreeVector p1v = E1 * dir1;

    // Angle between photon 1 and 2 using cosine law
    G4double cos_theta12 = (E3 * E3 - E1 * E1 - E2 * E2) / (2.0 * E1 * E2);
    if (std::abs(cos_theta12) > 1.0) {
        cos_theta12 = (cos_theta12 > 1.0) ? 1.0 : -1.0;
    }
    G4double theta12 = std::acos(cos_theta12);
    G4double sin_theta12 = std::sin(theta12);
    G4double phi2 = 2.0 * M_PI * G4UniformRand();

    // Direction for photon 2
    G4ThreeVector dir2(sin_theta12 * std::cos(phi2), sin_theta12 * std::sin(phi2), cos_theta12);
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