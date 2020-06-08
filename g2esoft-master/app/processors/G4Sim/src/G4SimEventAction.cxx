#include "G4SimEventAction.h"
#include "G4SimSteppingAction.h"

#include "TreeProc/LogStream.h"

#include "G4Event.hh"

using namespace TreeProc;
using namespace std;

G4SimEventAction::G4SimEventAction () {
}

G4SimEventAction::~G4SimEventAction () {
}

void G4SimEventAction::BeginOfEventAction (const G4Event *event) {
  log("debug") << "\x1b[33m[G4Sim] G4SimEventAction::BeginOfEventAction starts.\x1b[0m" << endl;

  if(_stepAction){
    _stepAction->ClearEvent();
    _stepAction->SetEventID(event->GetEventID());
  }

  log("debug") << "\x1b[33m[G4Sim] G4SimEventAction::BeginOfEventAction finished.\x1b[0m" << endl;
}

void G4SimEventAction::EndOfEventAction (const G4Event */*event*/) {
  log("debug") << "\x1b[33m[G4Sim] G4SimEventAction::EndOfEventAction starts.\x1b[0m" << endl;

  log("debug") << "\x1b[33m[G4Sim] G4SimEventAction::EndOfEventAction finished.\x1b[0m" << endl;
}
