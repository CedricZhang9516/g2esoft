#ifndef OBJECTS_H
#define OBJECTS_H

#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMatrixDSym.h>

#include <TRef.h>
#include <TRefArray.h>

#include "ref_vector.h"

namespace g2esoft{

  class MCParticle; // forward declaration
  
  class MCStep : public TObject{
  public:
    MCStep();
    MCStep(const MCStep &mcstep);
    ~MCStep(){}
    
  public:
    TVector3 _pos;
    double _time;
    TLorentzVector _p;
    TRef _mcp;

    // type-safe accessors for MCStep
    const MCParticle *getMCP()const;
    
    ClassDef(MCStep,4);
  };
    
    
  class MCParticle : public TObject{
  public:
    MCParticle();
    MCParticle(const MCParticle &mcp);
    ~MCParticle(){}
    
  public:
    int _pdg;
    TLorentzVector _p;
    TVector3 _pol;
    double _time;
    TVector3 _prodVertex;
    int _muonID;
    
    std::vector<MCParticle *> _daughters; //->
    MCParticle *_parent;
    
    const_ref_vector<MCStep> _steps;
    
    ClassDef(MCParticle,4);
  };

  class SimHit : public TObject{
  public:
    SimHit();
    SimHit(const SimHit &simhit);
    ~SimHit(){}

  public:
    TVector3 _pos[2];
    double _time;
    double _edep;
    unsigned int _type;

    TRef _mcStep;

    // type-safe accessors for MCStep
    const MCStep * getMCStep()const;

    ClassDef(SimHit,4);
  };

  class StripHit : public TObject{
  public:
    StripHit();
    StripHit(const StripHit &striphit);
    ~StripHit(){}

  public:
    unsigned int _stripID;
    unsigned int _tot;  // in nanosec
    double _time; // in nanosec
    double _edep;

    const StripHit *_prevHit;
    const StripHit *_nextHit;

    //const_ref_vector<SimHit> _simHits;
    std::vector<const SimHit *> _simHits;

    ClassDef(StripHit,6);
  };

  class StripCluster : public TObject{
  public:
    StripCluster();
    StripCluster(const StripCluster &stripclust);
    ~StripCluster(){}

    unsigned int vaneID() const{ return _stripID/0x1000000; }
    unsigned int sensorID() const{ return (_stripID%0x1000000)/0x10000; }
    bool isZStrip() const{ return (sensorID() & 0x04); }
    bool isRStrip() const{ return !isZStrip(); }

  public:
    unsigned int _stripID;
    double _pos;
    double _posErr;
    double _time;
    double _edep;

    //const_ref_vector<StripHit> _stripHits;
    std::vector<const StripHit *> _stripHits;

    ClassDef(StripCluster,5);
  };

  class RecoHit : public TObject{
  public:
    RecoHit();
    RecoHit(const RecoHit &recohit);
    ~RecoHit(){}

  public:
    TVector3 _pos;
    double _time;
    unsigned int _trueType;

    //const_ref_vector<StripCluster> _stripClusters;
    std::vector<const StripCluster *> _stripClusters;

    ClassDef(RecoHit,6);
  };

  enum RecoHitTrueTypeBit{
    IsGhostHit=0,
    IsNoiseHit=1,
  };

  class Track : public TObject{
  public:
    Track();
    Track(const Track &track);
    ~Track(){}

  public:
    TVector3 _pos;
    TVector3 _p;
    double _time;
    int _charge;
    double _chi2ndf;
    TMatrixDSym _covM;
    unsigned int _type;

    std::vector<const RecoHit *> _recoHits;
    TRef _mcp;
    TRef _mcStep;
    const MCParticle * getMCP()const;
    const MCStep *getMCStep()const;
    const Track *_prev;
    const Track *_next;

    ClassDef(Track,8);
  };

  enum TrackTypeBit{
    IsGhostTrack=0,
  };

}

#endif
