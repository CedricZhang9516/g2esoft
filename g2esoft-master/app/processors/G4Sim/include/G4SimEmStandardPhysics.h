#ifndef G4SimEmStandardPhysics_h
#define G4SimEmStandardPhysics_h 1

#include "G4VPhysicsConstructor.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4SimEmStandardPhysics : public G4VPhysicsConstructor
{
  public:
    G4SimEmStandardPhysics(const G4String& name = "standard");
   ~G4SimEmStandardPhysics();

  public:
    // This method is dummy for physics
    void ConstructParticle() {};

    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
    void ConstructProcess();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
