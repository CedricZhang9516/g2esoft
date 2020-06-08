// TrkFitting.h

#ifndef TrkFitting_h
#define TrkFitting_h

#include "TreeProc/Processor.h"
#include "app/objects.h"

#include <vector>
#include <string>

namespace genfit{
  class AbsKalmanFitter;
  class MeasurementCreator;
  class MaterialEffects;
  class EventDisplay;
  class Track;
  class DetPlane;
}

namespace g2esoft{
  class TrkFitting : public TreeProc::Processor {
  public:
    TrkFitting(const char *procName, const char *instName);
    virtual ~TrkFitting();
    
    virtual void config(TreeProc::Parameter *);
    virtual void process(int nEvent, TreeProc::Parameter *);
    virtual void finish(TreeProc::Parameter *);

    void setRecoHitPosition(genfit::Track *fitTrack, const RecoHit *recoHit);
    void setStripClusterPosition(genfit::Track *fitTrack, const StripCluster *stripCluster);
    genfit::DetPlane* getRecoHitPlane(const RecoHit *recoHit);
    genfit::DetPlane* getStripClusterPlane(const StripCluster *stripCluster);
    double getTrueTime(const RecoHit *recoHit, const bool doFirst=true);
    double getTrueTime(const StripCluster *stripCluster, const bool doFirst=true);

  private:
    std::vector<Track *> _tempTracks;
    std::vector<const Track*> _fittedTracks;

    std::string _parInputTrackName;
    std::string _parFittedTrackName;
    std::string _parFitterName;
    std::string _parDetectorGeometry;
    bool _parShowEventDisplay;
    bool _parUseStripCluster;
    bool _parRemoveGhostTrack;
    bool _parConnectTrack;
    double _parMomentumMin;
    double _parChi2NDFMax;
    double _parGhostTimeDiff;
    double _parConnectChi2NDFMax;
    double _parConnectTimeDiff;
    double _parConnectSearchTime;
    int _parNMaxIter;
    int _parNMinIter;
    int _parNMaxFailed;
    double _parDPVal;
    double _parDRelChi2;
    int _parNVane;
    std::vector<double> _parSensorOriginX;
    std::vector<double> _parSensorOriginZ;
    double _parPosFirstStrip;
    double _parStripPitch;
    double _parSensorDistance;
    double _parSensorSize;
    int _parNStripPerRow;
    int _parTimeDiv;

    std::string _parMagnetData;
    bool _parUniformField;
    bool _parInterpolateField;

    genfit::AbsKalmanFitter *_fitter;
    genfit::MeasurementCreator *_measurementCreator;
    genfit::MaterialEffects *_matEff;
    genfit::EventDisplay *_display;

  };
}

#endif
