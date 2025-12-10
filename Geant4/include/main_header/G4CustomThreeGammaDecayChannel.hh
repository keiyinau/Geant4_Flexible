// G4CustomThreeGammaDecayChannel.hh
#ifndef G4CUSTOMTHREEGAMMADECAYCHANNEL_HH
#define G4CUSTOMTHREEGAMMADECAYCHANNEL_HH

#include "G4VDecayChannel.hh"

class G4CustomThreeGammaDecayChannel : public G4VDecayChannel {
public:
    G4CustomThreeGammaDecayChannel(const G4String& theParentName, G4double theBR, G4bool useLoops = false);
    virtual ~G4CustomThreeGammaDecayChannel() {}

    virtual G4DecayProducts* DecayIt(G4double parentMass);

private:
    G4bool loops;
};

#endif