#include "app/objects.h"

//ClassImp(g2esoft::MCParticle)
//ClassImp(g2esoft::SimHit)

namespace g2esoft{

  MCStep::MCStep():
    _time(0.)
  {
  }

  MCStep::MCStep(const MCStep &mcstep):
    _pos (mcstep._pos),
    _time(mcstep._time),
    _p   (mcstep._p)
  {
  }

  const MCParticle *MCStep::getMCP()const
  {
    if(_mcp.IsValid()){
      return static_cast<MCParticle *>(_mcp.GetObject());
    }else{
      return 0;
    }
  }

  MCParticle::MCParticle():
    _pdg   (0),
    _time  (0.),
    _muonID(0),
    _parent(nullptr)
  {
  }

  MCParticle::MCParticle(const MCParticle &mcp):
    _pdg       (mcp._pdg),
    _p         (mcp._p),
    _pol       (mcp._pol),
    _time      (mcp._time),
    _prodVertex(mcp._prodVertex),
    _muonID    (mcp._muonID),
    _parent    (nullptr)
  {
  }

  SimHit::SimHit():
    _time(0.),
    _edep(0.),
    _type(0)
  {
  }

  SimHit::SimHit(const SimHit &simhit):
    _pos {simhit._pos[0],simhit._pos[1]},
    _time(simhit._time),
    _edep(simhit._edep),
    _type(simhit._type)
  {
  }

  const MCStep * SimHit::getMCStep()const
  {
    if(_mcStep.IsValid()){
      return static_cast<MCStep *>(_mcStep.GetObject());
    }else{
      return 0;
    }
  }

  StripHit::StripHit():
    _stripID(0),
    _tot(0),
    _time(0.),
    _prevHit(nullptr),
    _nextHit(nullptr)
  {
  }

  StripHit::StripHit(const StripHit &striphit):
    _stripID(striphit._stripID),
    _tot(striphit._tot),
    _time(striphit._time),
    _prevHit(nullptr),
    _nextHit(nullptr)
  {
  }

  StripCluster::StripCluster():
    _stripID(0),
    _pos(0.),
    _posErr(0.),
    _time(0.),
    _edep(0.)
  {
  }

  StripCluster::StripCluster(const StripCluster &stripclust):
    _stripID(stripclust._stripID),
    _pos(stripclust._pos),
    _posErr(stripclust._posErr),
    _time(stripclust._time),
    _edep(stripclust._edep)
  {
  }
    

  RecoHit::RecoHit():
    _time(0.),
    _trueType(0)
  {
  }

  RecoHit::RecoHit(const RecoHit &recohit):
    _pos(recohit._pos),
    _time(recohit._time),
    _trueType(recohit._trueType)
  {
  }

  Track::Track():
    _time(0.),
    _charge(0),
    _chi2ndf(0),
    _type(0),
    _prev(nullptr),
    _next(nullptr)
  {
  }

  Track::Track(const Track &track):
    _pos(track._pos),
    _p(track._p),
    _time(track._time),
    _charge(track._charge),
    _chi2ndf(0),
    _type(track._type),
    _prev(nullptr),
    _next(nullptr)
  {
  }

  const MCParticle *Track::getMCP()const
  {
    if(_mcp.IsValid()){
      return static_cast<MCParticle *>(_mcp.GetObject());
    }else{
      return 0;
    }
  }

  const MCStep *Track::getMCStep()const
  {
    if(_mcStep.IsValid()){
      return static_cast<MCStep *>(_mcStep.GetObject());
    }else{
      return 0;
    }
  }

}

