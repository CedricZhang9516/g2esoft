// G4Sim.h

#ifndef G4Sim_h
#define G4Sim_h 1

#include "TreeProc/Processor.h"
#include "app/objects.h"

class G4SimRunManager;

namespace g2esoft{
  class G4Sim : public TreeProc::Processor {
  public:
    G4Sim (const char *procName, const char *instName);
    virtual ~G4Sim ();

    virtual void config (TreeProc::Parameter *);
    virtual void process (int nEvent, TreeProc::Parameter *);
    virtual void finish (TreeProc::Parameter *);

  private:
    int _parNEvents;

    std::vector<MCParticle *> _mcParticles;
    std::vector<const MCParticle *> _mcParticlesBuffer;
    std::vector<const MCStep *> _mcSteps;
    std::vector<const SimHit *> _simHits;

    G4SimRunManager *_runManager;
  };
}

#endif
