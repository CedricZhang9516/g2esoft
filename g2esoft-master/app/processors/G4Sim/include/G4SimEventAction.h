#ifndef G4SIMEVENTACTION_H
#define G4SIMEVENTACTION_H

#include "G4UserEventAction.hh"

class G4Event;
class G4SimSteppingAction;

class G4SimEventAction : public G4UserEventAction {
  public:
    G4SimEventAction ();
    virtual ~G4SimEventAction ();

    virtual void BeginOfEventAction (const G4Event *);
    virtual void EndOfEventAction (const G4Event *);

    void SetSteppingAction(G4SimSteppingAction *stepAction){_stepAction=stepAction;}

  private:
    G4SimSteppingAction *_stepAction;
};

#endif
