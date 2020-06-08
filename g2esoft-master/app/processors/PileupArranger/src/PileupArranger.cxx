// PileupArranger.cxx

#include "PileupArranger.h"
#include "TreeProc/TreeManager.h"
//#include "g2track/g2track.h"
#include "TreeProc/Factory.h"
#include "TreeProc/RunManager.h"

#include <iostream>
#include <unordered_map>
#include "TRandom3.h"

using namespace g2esoft;
using namespace TreeProc;
using namespace std;

static Factory<EventArranger, PileupArranger> _f("PileupArranger");

PileupArranger::PileupArranger(const char *procName, const char *instName) 
  : EventArranger(procName, instName),
    _random3(0)
{
}

PileupArranger::~PileupArranger()
{
}

void PileupArranger::config(Parameter * param){
  _parSimHitName = param->get("SimHitName",string("PileupSimHits"));
  TreeManager::instance()->registerBranch(_parSimHitName.c_str(), &_simHits);
  _parMCParticleName = param->get("MCParticleName",string("PileupMCParticles"));
  TreeManager::instance()->registerBranch(_parMCParticleName.c_str(), &_mcps);
  _parMCStepName = param->get("MCStepName",string("PileupMCSteps"));
  TreeManager::instance()->registerBranch(_parMCStepName.c_str(), &_mcSteps);

  _parInputSimHitName = param->get("InputSimHitName",string("SimHits"));
  _parInputMCParticleName = param->get("InputMCParticleName",string("MCParticles"));
  _parInputMCStepName = param->get("InputMCStepName",string("MCSteps"));

  _parModifyTime = param->get("ModifyTime", (bool)true);
  _parTimeRange = param->get("TimeRange",(double)10.0);
  _parNumberToCombine = param->get("NumberToCombine",10);
  _parRandomSeed = param->get("RandomSeed",(int)0);

  if(_parModifyTime){
    _random3 = new TRandom3((unsigned int)_parRandomSeed);
  }
}

bool PileupArranger::arrange(int nEventIn, Parameter *)
{
  for(unsigned int i=0; i<_mcps.size(); i++){
    delete _mcps[i];
  }
  _mcps.clear();
  for(unsigned int i=0; i<_simHits.size(); i++){
    delete _simHits[i];
  }
  _simHits.clear();
  for(unsigned int i=0; i<_mcSteps.size(); i++){
    delete _mcSteps[i];
  }
  _mcSteps.clear();
  
  std::unordered_map<const MCParticle*, MCParticle*> mcpMap;
  std::unordered_map<const MCStep*, MCStep*> mcstepMap;

  for(unsigned int i=0;i<_parNumberToCombine;i++){
    // loading the next entry
    int nev = RunManager::instance()->getNextEntry();
        
    // determine new time
    double newTime = 0.;
    if(_parModifyTime) _random3->Rndm()*_parTimeRange;
    double thisTime = 0.;
    
    // MCStep
    const vector<const MCStep *> &mcsteps = TreeManager::instance()->getBranchVec<MCStep>(_parInputMCStepName.c_str());
    for(unsigned int nmcst = 0;nmcst < mcsteps.size(); nmcst ++){
      MCStep *mcstep = new MCStep(*(mcsteps[nmcst]));
      if(nmcst==0 && _parModifyTime) thisTime = mcstep->_time;
      mcstep->_time += (newTime-thisTime);
      _mcSteps.push_back(mcstep);
      mcstepMap[mcsteps[nmcst]] = mcstep;
    }

    // MCParticle
    const vector<const MCParticle *> &mcps = TreeManager::instance()->getBranchVec<MCParticle>(_parInputMCParticleName.c_str());
    for(unsigned int nmc = 0;nmc < mcps.size(); nmc ++){
      MCParticle *mcp = new MCParticle(*(mcps[nmc]));
      mcp->_time += (newTime-thisTime);
      _mcps.push_back(mcp);
      mcpMap[mcps[nmc]] = mcp;
    }

    // link parent, daughters, MCParticle and MCStep
    for(unsigned int nmc = 0;nmc < mcps.size(); nmc ++){
      if( mcps[nmc]->_parent ){
	mcpMap.at(mcps[nmc])->_parent = mcpMap.at(mcps[nmc]->_parent);
      }
      for(unsigned int ndaughter = 0; ndaughter < mcps[nmc]->_daughters.size(); ndaughter++ ){
	mcpMap.at(mcps[nmc])->_daughters.push_back( mcpMap.at(mcps[nmc]->_daughters.at(ndaughter)) );
      }
      for(unsigned int nmcstep = 0; nmcstep<mcps[nmc]->_steps.size(); nmcstep++ ){
	mcpMap.at(mcps[nmc])->_steps.push_back( mcstepMap.at(mcps[nmc]->_steps.at(nmcstep)) );
	mcstepMap.at(mcps[nmc]->_steps.at(nmcstep))->_mcp = mcpMap.at(mcps[nmc]);
      }
    }

    // SimHit
    const vector<const SimHit *> &simhits = TreeManager::instance()->getBranchVec<SimHit>(_parInputSimHitName.c_str());
    for(unsigned int nsim = 0;nsim < simhits.size(); nsim ++){
      SimHit *outSimHit = new SimHit(*(simhits[nsim]));
      outSimHit->_time += (newTime-thisTime);
      outSimHit->_mcStep = mcstepMap.at(simhits[nsim]->getMCStep());
      
      _simHits.push_back(outSimHit);
    }
    
  }
  return true;  
}
