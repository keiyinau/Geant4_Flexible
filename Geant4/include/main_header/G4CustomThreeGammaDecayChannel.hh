#ifndef G4CUSTOM_THREE_GAMMA_DECAY_CHANNEL_HH
#define G4CUSTOM_THREE_GAMMA_DECAY_CHANNEL_HH

#include "G4VDecayChannel.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

class G4CustomThreeGammaDecayChannel : public G4VDecayChannel {
public:
    G4CustomThreeGammaDecayChannel(const G4String& theParentName, G4double theBR);
    virtual ~G4CustomThreeGammaDecayChannel() {}
    virtual G4DecayProducts* DecayIt(G4double parentMass = 0.0);
};

#endif