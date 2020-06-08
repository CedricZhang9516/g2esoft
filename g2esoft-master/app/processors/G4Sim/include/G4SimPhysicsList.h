#ifndef G4SimPhysicsList_h
#define G4SimPhysicsList_h 1

#include <vector>

#include "G4VModularPhysicsList.hh"

class G4VPhysicsConstructor;
class G4SimStepMax;

class G4SimPhysicsList: public G4VModularPhysicsList {
public:
  G4SimPhysicsList ();
  virtual ~G4SimPhysicsList ();

  void ConstructParticle();
  void ConstructProcess ();

  void SetCuts ();
  void SetCutForGamma (G4double);
  void SetCutForElectron (G4double);
  void SetCutForPositron (G4double);

  void SetStepMax (G4double);
  void AddStepMax ();

  /// Make sure that the EM physics list is empty.
  void ClearEMPhysics ();

  /// Make sure that the hadron physics list is empty.
  void ClearHadronPhysics ();


private:
  typedef std::vector<G4VPhysicsConstructor *> PhysicsListVector;

  G4double _cutForGamma;
  G4double _cutForElectron;
  G4double _cutForPositron;

  G4VPhysicsConstructor *_decayPhysicsList;
  G4VPhysicsConstructor *_emPhysicsList;

  PhysicsListVector *_emPhysicsListVector;
  PhysicsListVector *_hadronPhysicsListVector;

  G4SimStepMax *_stepMaxProcess;
};
#endif
