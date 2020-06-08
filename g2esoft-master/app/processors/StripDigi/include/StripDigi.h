// StripDigi.h

#ifndef StripDigi_h
#define StripDigi_h

#include "TreeProc/Processor.h"
#include "app/objects.h"

#include <vector>
#include <string>

class TSpline3;

namespace g2esoft{
  class StripDigi : public TreeProc::Processor {
  public:
    StripDigi(const char *procName, const char *instName);
    virtual ~StripDigi();
    
    virtual void config(TreeProc::Parameter *);
    virtual void process(int nEvent, TreeProc::Parameter *);
    virtual void finish(TreeProc::Parameter *);
    
  private:
    std::vector<const StripHit *> _stripHits;

    bool _parStripInG4;
    std::string _parStripHitName;
    std::string _parSimHitName;
    int  _parNVane;
    int  _parNStripPerRow;
    bool _parDoubleRow;
    std::vector<double> _parSensorOriginX;
    std::vector<double> _parSensorOriginZ;
    double _parSensorSize;
    double _parPosFirstStrip;
    double _parStripPitch;
    double _parInactiveEdge;
    double _parInactiveMiddle;
    int    _parTimeDiv;
    double _parEdepMIPThreshold;
    double _parEdepMIPThreshold2;
    std::string _parTimeDigiMethod;
    std::string _parCRRCParamFileName;
    std::string _parDiffParamFileName;

    std::vector<TSpline3 *> _spCRRCParam;
    std::vector<TSpline3 *> _spDiffParam;
  };
}

#endif
