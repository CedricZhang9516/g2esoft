#ifndef G4SIMSTEPPINGACTION_H
#define G4SIMSTEPPINGACTION_H 1

#include <unordered_map>

#include "app/objects.h"

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"

class G4SimSteppingAction : public G4UserSteppingAction {
  public:
  G4SimSteppingAction (std::vector<g2esoft::MCParticle *> &mcParticles, std::vector<const g2esoft::MCStep *> &mcSteps, std::vector<const g2esoft::SimHit *> &simHits);
    virtual ~G4SimSteppingAction ();

    virtual void UserSteppingAction (const G4Step *);
    void ClearEvent();
    void SetEventID(const G4int eventID){_currentEventID = eventID;}
    void CreateHit();

  private:
    std::vector<g2esoft::MCParticle *> &_mcParticles;
    std::vector<const g2esoft::MCStep *> &_mcSteps;
    std::vector<const g2esoft::SimHit *> &_simHits;
    std::unordered_map<G4int, g2esoft::MCParticle *> _mcpMap;
 
    bool _saveStep;

    G4int _currentEventID;
    G4int _currentTrackID;
    G4int _currentPdg;
    G4double _edepSum;
    G4ThreeVector _prePosition;
    G4ThreeVector _postPosition;
    G4ThreeVector _preMomentum;
    G4ThreeVector _prePolarization;
    G4double _preEnergy;
    G4double _preTime;
};

#endif
