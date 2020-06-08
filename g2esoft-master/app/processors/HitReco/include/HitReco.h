// HitReco.h

#ifndef HitReco_h
#define HitReco_h

#include "TreeProc/Processor.h"
#include "app/objects.h"

#include <vector>
#include <string>

namespace g2esoft{
  class HitReco : public TreeProc::Processor {
  public:
    HitReco(const char *procName, const char *instName);
    virtual ~HitReco();
    
    virtual void config(TreeProc::Parameter *);
    virtual void process(int nEvent, TreeProc::Parameter *);
    virtual void finish(TreeProc::Parameter *);
    
  public:
    TVector3 stripPosition(const unsigned int stripID);
    TVector3 stripPosition(const unsigned int stripID, const double stripPos);
    inline bool isZStrip(const unsigned int stripID);
    inline bool isRStrip(const unsigned int stripID);
    inline unsigned int vaneID(const unsigned int stripID);
    inline unsigned int sensorID(const unsigned int stripID);
    inline unsigned int stripOnlyID(const unsigned int stripID);

  private:
    std::vector<const RecoHit *> _recoHits;
    
    std::string _parRecoHitName;
    std::string _parStripClusterName;

    // parametes for hit reconstruction
    int    _parNVane;
    bool   _parDoubleRow;
    std::vector<double> _parSensorOriginX;
    std::vector<double> _parSensorOriginZ;
    double _parStripPitch;
    double _parSensorSize;
    double _parInactiveMiddle;

    int _parTimeDiv;
    double _parTimeMatch;

  };
}

#endif
