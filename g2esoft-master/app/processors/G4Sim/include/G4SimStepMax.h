#ifndef G4SimStepMax_h
#define G4SimStepMax_h 1

#include "G4VDiscreteProcess.hh"

class G4ParticleDefinition;
class G4Step;

class G4SimStepMax : public G4VDiscreteProcess
{
  public:
    G4SimStepMax (const G4String &processName = "UserStepMax");
    G4SimStepMax (G4SimStepMax &);

    ~G4SimStepMax ();

    G4bool IsApplicable (const G4ParticleDefinition &);

    void SetStepMax (G4double);

    G4double GetStepMax () { return _maxChargedStep; };

    G4double PostStepGetPhysicalInteractionLength (const G4Track &track, G4double previousStepSize, G4ForceCondition *condition);

    G4VParticleChange *PostStepDoIt (const G4Track &, const G4Step &);

  protected:
    G4double GetMeanFreePath (const G4Track &, G4double, G4ForceCondition *);

  private:
    // hide assignment operator as private
    G4SimStepMax & operator= (const G4SimStepMax &right);
    G4SimStepMax (const G4SimStepMax &);

  private:
    G4double _maxChargedStep;

};

#endif
