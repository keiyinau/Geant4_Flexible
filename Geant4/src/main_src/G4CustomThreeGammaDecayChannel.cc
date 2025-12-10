// G4CustomThreeGammaDecayChannel.cc
#include "G4CustomThreeGammaDecayChannel.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4LorentzVector.hh"
#include "G4Gamma.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>

G4CustomThreeGammaDecayChannel::G4CustomThreeGammaDecayChannel(const G4String& theParentName, G4double theBR, G4bool useLoops)
    : G4VDecayChannel("Custom 3 Gamma Decay", theParentName, theBR, 3, "gamma", "gamma", "gamma"), loops(useLoops) {}

G4DecayProducts* G4CustomThreeGammaDecayChannel::DecayIt(G4double parentMass) {
    if (parentMass <= 0.0) parentMass = GetParentMass();

    // Get o-Ps polarization from dynamic particle (parent is static definition; use track pol if available, or assume from setup)
    // Note: In decay, polarization comes from the track; here assume S from previous propagation (e.g., set to [0,0,1] for demo if not available)
    G4ThreeVector S(0., 0., 1.);  // Placeholder; in full sim, pass from G4Track via custom hook or assume longitudinal

    G4double M = parentMass;  // o-Ps mass ~1.022 MeV
    const static G4double min_energy = 0.001 * M / 2;  // IR cutoff
    const G4int n_points = 1;  // Generate one event per call
    const G4int n_probe = 1000000;  // Number of probes to estimate max density

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Updated dalitz_density with loops switch
    auto dalitz_density = [M, S, this](G4double e1, G4double e2, G4double e3, G4ThreeVector k1, G4ThreeVector k2, G4ThreeVector k3) -> G4double {
        G4double z = M - e1 - e2;
        if (e1 < min_energy || e2 < min_energy || z < min_energy || e1 >= M/2 || e2 >= M/2 || z >= M/2) return 0.0;

        G4double term1 = std::pow((M/2 - e1), 2) / (std::pow(e2, 2) * std::pow(z, 2));
        G4double term2 = std::pow((M/2 - e2), 2) / (std::pow(e1, 2) * std::pow(z, 2));
        G4double term3 = std::pow((M/2 - z), 2) / (std::pow(e1, 2) * std::pow(e2, 2));
        G4double unpol = term1 + term2 + term3;  // Ore-Powell tree-level unpol

        G4double density;
        if (!loops) {  // Tree-level: Unpol + simple pol corr
            G4ThreeVector cross = k1.cross(k2);
            G4double asymmetry=0.;
            G4double polTerm = S.dot(cross.unit()) *asymmetry;  // Asymmetry ~5% (tunable; from Phys. Rev. A 81, 2010)
            density = unpol * (1.0 + polTerm);  // Tree-level polarized rate
        } else {  // Loop-level: One-loop amplitudes A1,A2,A3 (from Phys. Rev. A 72, 032501, 2005)
            // Analytic A1,A2,A3 (simplified; full from hep-ph/0506213)
            G4double s12 = 2 * (e1 * e2 - k1.dot(k2));  // Mandelstam-like (2p1·p2)
            G4double s13 = 2 * (e1 * e3 - k1.dot(k3));
            G4double s23 = 2 * (e2 * e3 - k2.dot(k3));
            G4double m2 = electron_mass_c2 * electron_mass_c2;

            // One-loop corrected A_i (approx from arXiv:hep-ph/0506213; exact eval needs full code)
            G4double A1 = (1.0 / (s12 * s13)) * (log(s12 / m2) - log(s13 / m2)) + fine_structure_const / pi;  // O(α) example
            G4double A2 = (1.0 / (s23 * s13)) * (log(s23 / m2) - log(s13 / m2)) + fine_structure_const / pi;
            G4double A3 = (1.0 / (s12 * s23)) * (log(s12 / m2) - log(s23 / m2)) + fine_structure_const / pi;

            // Full |M|^2 from tensor (simplified; per Phys. Rev. Lett. 76, 4903, 1996)
            G4double M_sq = 2 * (A1 * A1 * s23 * s23 + A2 * A2 * s13 * s13 + A3 * A3 * s12 * s12) 
                          - (A1 * A2 * s12 * s13 * s23 / m2);  // Bose-symmetric terms

            // Add pol corr at loop level (from hep-ph/9708450)
            G4ThreeVector cross = k1.cross(k2);
            G4double polTerm = S.dot(cross.unit()) * (fine_structure_const / pi) * 0.01;  // Loop-suppressed asymmetry
            density = M_sq * (1.0 + polTerm);  // Loop-corrected polarized
        }
        return density;
    };

    // Find approximate max density for rejection sampling
    G4double max_d = 0.0;
    for (G4int i = 0; i < n_probe; ++i) {
        G4double x = min_energy + dis(gen) * (M / 2 - 2 * min_energy);
        G4double y = min_energy + dis(gen) * (M / 2 - 2 * min_energy);
        G4double z = M - x - y;
        G4ThreeVector dummy_k1(1,0,0), dummy_k2(0,1,0), dummy_k3(0,0,1);  // Dummy dirs for max est
        if (z > min_energy && z < M / 2 && x + y < M) {
            G4double d = dalitz_density(x, y, z, dummy_k1, dummy_k2, dummy_k3);
            if (d > max_d) max_d = d;
        }
    }
    max_d *= 1.1;  // Safety margin

    // Generate one event using rejection sampling
    G4double E1, E2, E3;
    G4ThreeVector dir1, dir2, dir3;
    G4int accepted = 0;
    while (accepted == 0) {
        G4double x = min_energy + dis(gen) * (M / 2 - min_energy);
        G4double y = min_energy + dis(gen) * (M / 2 - min_energy);
        E3 = M - x - y;
        if (E3 > min_energy && E3 < M / 2 && x + y < M) {
            // Sample directions
            G4double cost = 2.0 * dis(gen) - 1.0;
            G4double sint = std::sqrt(1.0 - cost * cost);
            G4double phi = 2.0 * M_PI * dis(gen);
            dir1 = G4ThreeVector(sint * std::cos(phi), sint * std::sin(phi), cost);
            G4ThreeVector k1 = x * dir1;

            G4double cos_theta12 = (E3 * E3 - x * x - y * y) / (2.0 * x * y);
            if (std::abs(cos_theta12) > 1.0) cos_theta12 = (cos_theta12 > 1.0) ? 1.0 : -1.0;
            G4double sin_theta12 = std::sqrt(1.0 - cos_theta12 * cos_theta12);
            G4double phi2 = 2.0 * M_PI * dis(gen);

            dir2 = G4ThreeVector(sin_theta12 * std::cos(phi2), sin_theta12 * std::sin(phi2), cos_theta12);
            dir2.rotateUz(dir1);
            G4ThreeVector k2 = y * dir2;

            G4ThreeVector k3 = -k1 - k2;
            dir3 = k3.unit();

            G4double d = dalitz_density(x, y, E3, k1, k2, k3);
            if (dis(gen) * max_d < d) {
                // Shuffle energies/dirs to avoid bias
                std::vector<G4double> E_shuffle = {x, y, E3};
                std::vector<G4ThreeVector> dir_shuffle = {dir1, dir2, dir3};
                std::shuffle(E_shuffle.begin(), E_shuffle.end(), gen);
                std::shuffle(dir_shuffle.begin(), dir_shuffle.end(), gen);
                E1 = E_shuffle[0]; E2 = E_shuffle[1]; E3 = E_shuffle[2];
                dir1 = dir_shuffle[0]; dir2 = dir_shuffle[1]; dir3 = dir_shuffle[2];
                accepted = 1;
            }
        }
    }

    // Create parent at rest
    G4DynamicParticle* parentParticle = new G4DynamicParticle(GetParent(), G4ThreeVector(0., 0., 0.), 0.0);
    G4DecayProducts* products = new G4DecayProducts(*parentParticle);
    delete parentParticle;

    // Create photons with directions and polarizations
    G4DynamicParticle* gamma1 = new G4DynamicParticle(G4Gamma::Gamma(), E1 * dir1);
    G4DynamicParticle* gamma2 = new G4DynamicParticle(G4Gamma::Gamma(), E2 * dir2);
    G4DynamicParticle* gamma3 = new G4DynamicParticle(G4Gamma::Gamma(), E3 * dir3);

    // Set polarizations (linear, perp to dir; random phi)
    G4double phi_pol = twopi * dis(gen);
    G4ThreeVector pol1(std::cos(phi_pol), std::sin(phi_pol), 0.);
    pol1.rotateUz(dir1);
    gamma1->SetPolarization(pol1);

    phi_pol = twopi * dis(gen);
    G4ThreeVector pol2(std::cos(phi_pol), std::sin(phi_pol), 0.);
    pol2.rotateUz(dir2);
    gamma2->SetPolarization(pol2);

    phi_pol = twopi * dis(gen);
    G4ThreeVector pol3(std::cos(phi_pol), std::sin(phi_pol), 0.);
    pol3.rotateUz(dir3);
    gamma3->SetPolarization(pol3);

    products->PushProducts(gamma1);
    products->PushProducts(gamma2);
    products->PushProducts(gamma3);

#ifdef G4VERBOSE
    if (GetVerboseLevel() > 1) {
        G4cout << "G4CustomThreeGammaDecayChannel::DecayIt - Created 3 photons with energies: "
               << E1 / keV << ", " << E2 / keV << ", " << E3 / keV << " keV (loops=" << loops << ")" << G4endl;
    }
#endif

    return products;
}