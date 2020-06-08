#ifndef G2ESOFT_PILEUPARRANGER_H
#define G2ESOFT_PILEUPARRANGER_H

#include "TreeProc/EventArranger.h"
#include "app/objects.h"

#include <vector>
#include <string>

class TRandom3;

namespace g2esoft{
  class PileupArranger : public TreeProc::EventArranger {
  public:
    PileupArranger(const char *procName, const char *instName);
    virtual ~PileupArranger();
    
    virtual void config(TreeProc::Parameter *);
    virtual bool arrange(int nEventIn, TreeProc::Parameter *);
    virtual void finish(TreeProc::Parameter *){}
  private:
    
    std::vector<const SimHit *> _simHits;
    std::vector<const MCParticle *> _mcps;
    std::vector<const MCStep *> _mcSteps;
    
    std::string _parSimHitName;
    std::string _parMCParticleName;
    std::string _parMCStepName;

    std::string _parInputSimHitName;
    std::string _parInputMCParticleName;
    std::string _parInputMCStepName;

    bool _parModifyTime; // whether modify time or not
    double _parTimeRange; // time range of overlay
    int _parNumberToCombine; // number of pileup
    int _parRandomSeed;

    TRandom3 *_random3;
  };
}

#endif
//    TrackMultiplier(const int multiMode, const double range, const unsigned long seed=0);
//    void multiply(const vector<MCParticle *> *mcparts, const vector<SimHit *> *simHits, 
//		  vector<MCParticle *> *multiMCParts, vector<SimHit *> *multiSimHits);
