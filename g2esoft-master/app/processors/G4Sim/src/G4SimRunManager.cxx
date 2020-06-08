#include "G4SimRunManager.h"

void G4SimRunManager::BeamOnStart(G4int n_event, const char *macroFile, G4int n_select)
{
  _beamOnConfirmed = ConfirmBeamOnCondition();
  if(_beamOnConfirmed){
    numberOfEventToBeProcessed = n_event;
    ConstructScoringWorlds();
    RunInitialization();
    if(n_event>0){
      InitializeEventLoop(n_event,macroFile,n_select);
    }
  }
}

void G4SimRunManager::BeamOnNext(G4int i_event)
{
  if(_beamOnConfirmed && numberOfEventToBeProcessed>0){
    ProcessOneEvent(i_event);
    TerminateOneEvent();
    if(runAborted){
      BeamOnEnd();
      _beamOnConfirmed = false;
    }
  }
}

void G4SimRunManager::BeamOnEnd()
{
  TerminateEventLoop();
  RunTermination();
}
