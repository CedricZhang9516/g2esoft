#ifndef G4SIMRUNMANAGER_H
#define G4SIMRUNMANAGER_H

#include "G4RunManager.hh"

class G4SimRunManager : public G4RunManager {
 public:
  void BeamOnStart(G4int n_event, const char *macroFile=0, G4int n_select=1);
  void BeamOnNext(G4int i_event);
  void BeamOnEnd();

 protected:
  G4bool _beamOnConfirmed;
};

#endif
