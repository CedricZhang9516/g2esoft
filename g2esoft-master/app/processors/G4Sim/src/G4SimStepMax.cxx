#include "G4SimStepMax.h"

#include "TreeProc/LogStream.h"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"


using namespace TreeProc;
using namespace std;

G4SimStepMax::G4SimStepMax (const G4String &aName) : G4VDiscreteProcess (aName), _maxChargedStep (DBL_MAX) {
  if (verboseLevel > 0) {
   log("debug") << GetProcessName () << " is created "<< endl;
  }
}

G4SimStepMax::~G4SimStepMax () { }

G4SimStepMax::G4SimStepMax (G4SimStepMax &right) : G4VDiscreteProcess (right) { }

G4bool G4SimStepMax::IsApplicable (const G4ParticleDefinition &particle) {
  return (particle.GetPDGCharge () != 0.0);
}

void G4SimStepMax::SetStepMax (G4double step) { _maxChargedStep = step ; }

G4double G4SimStepMax::PostStepGetPhysicalInteractionLength (const G4Track &, G4double, G4ForceCondition *condition) {
  // condition is set to "Not Forced"
  *condition = NotForced;

  G4double ProposedStep = DBL_MAX;

  if (_maxChargedStep > 0.0)
    ProposedStep = _maxChargedStep ;

  return ProposedStep;
}

G4VParticleChange* G4SimStepMax::PostStepDoIt (const G4Track &aTrack, const G4Step &) {
  // do nothing
  aParticleChange.Initialize (aTrack);
  return &aParticleChange;
}

G4double G4SimStepMax::GetMeanFreePath (const G4Track &,G4double, G4ForceCondition *) {
  return 0.0;
}
