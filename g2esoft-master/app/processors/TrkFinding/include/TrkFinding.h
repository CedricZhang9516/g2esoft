// TrkFinding.h

#ifndef TrkFinding_h
#define TrkFinding_h

#include "TreeProc/Processor.h"
#include "app/objects.h"

#include <vector>
#include <unordered_map>

#include <TROOT.h>
#include <TF1.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TArc.h>
#include <TLine.h>
#include <TGraph.h>
#include <TMarker.h>
#include <TArrow.h>

namespace g2esoft{
  class TrkFinding : public TreeProc::Processor {
  public:
    TrkFinding(const char *name, const char *instName);
    virtual ~TrkFinding();
    
    virtual void config(TreeProc::Parameter *);
    virtual void process(int nEvent, TreeProc::Parameter *);
    virtual void finish(TreeProc::Parameter *);

  public:
    // function
    unsigned int getVaneID( double phi ); // from 0 to _geom_nvane-1. vane-ID = 0 at (x,y)=(+1,0). vane-ID is inclemented according to clockwise direction.
    double       calcPhi   ( double y, double x ); // from 0 to 2pi. phi=0 at (x,y)=(1,0). phi is increased according to anti-clockwise direction.
    void         inputHits( int index, const RecoHit* recohit );
    void         clearEvent(); // clear object at the end of each event
    void         clearTime (); // clear object at the end of each time-window
    void         printParameters();
    void         printHitInfo(int mode=0);
    void         printHitInfoVaneIDOrder     (){ printHitInfo(1); }
    void         printHitInfoRecTimeOrder    (){ printHitInfo(2); }
    void         printHitInfoMCTruthTimeOrder(){ printHitInfo(3); }
    void         printHitInfoZOrder          (){ printHitInfo(4); }
    void         calcOrder();
    void         houghTransformPhiZ();
    int          houghFitPhiZ();
    int          calcHoughResidualPhiZ();

    int   clusteringShort      ();
    int   clustering3DShort    ();
    int   clustering3DShortOld ();
    void  expandClustering     ( std::vector<int> &extrap_index, bool isforward );
    int   clusteringNextTimeBin( double time_window_min, double time_window_max );
    void  judgeFinding         ( std::vector<std::pair<int,int> > *result );
    void  evaluateTrkFinding   ();
    
    int   getNhits          (){ return _hitInfoRecoHits.size();     }
    int   getLineNoMax      (){ return _houghPhiZLineNoMax;         }
    int   getNLines         (){ return _houghPhiZFunc.size();       }
    int   getNActiveLines   (){ return _houghPhiZFunc.size();       }
    int   getClusterNoMax   (){ return _clusterNoMax;               }
    int   getNClusters      (){ return _mapClusteringResult.size(); }
    int   getNActiveClusters(){ return _mapActiveClusterNo.size();  }
    bool  isActiveClusterNo ( int clusterNo ){ return _mapActiveClusterNo.count(clusterNo); }
    
  private:
    const static int _cycleColor = 7;
    // canvas
    TCanvas* _canvas1;

    // geometry
    int     _parGeomNVane;
    double  _geomRInner;
    double  _geomROuter;
    double  _geomRPole;
    double  _geomZMin;
    double  _geomZMax;
    double* _geomPhi;

    // setting of hough-transformation and clustering
    int    _parCutHoughPhiZPeak;           // minimum entry to search hough-peak for  normal -track (used in hough peak search)
    int    _parCutHoughPhiZPeakDeltaRay;   // minimum entry to search hough-peak for deltaray-track (used in hough peak search)
    double _parCutTrackRangeThetaMin;      // range for  normal -track (used in hough peak search)
    double _parCutTrackRangeThetaMax;      // range for  normal -track (used in hough peak search)
    double _parCutDeltaRayRangeThetaMax;   // range for deltaray-track (used in hough peak search)

    double _parCutHoughPhiZSeedResi;         // tolerance range in cluster seed search
    double _parCutHoughPhiZSeedResiDeltaRay; // tolerance range in cluster seed search

    int    _parCutExtrapMiss;    // used in clustering
    int    _parCutExtrapNoCross; // used in clustering

    double _parCutExtrapTolerance;              // tolerance range in extrapolation
    double _parCutExtrapToleranceCoeffDphi;     // tolerance range in extrapolation : tolerance is changed as a function of delta phi

    //double _parEvtDispThresholdSignalEnergy;  // definition of signal energy in event display. It is related only to event display. (It does NOT contribute to the result of track finding.)
    double _parEvtDispThresholdMomentum;  // threshold of track momentum in event display. It is related only to event display. (It does NOT contribute to the result of track finding.)

    // hough transformation
    int    _parHoughNStepTheta;
    int    _parHoughNStepRho;
    double _parHoughRhoMin;
    double _parHoughRhoMax;

    int    _parThresholdSuccess;

    // parameter from xml
    // class name
    std::string _parRecoHitName;
    std::string _parTrackName;
    std::string _parMCParticleName;

    // flag
    int  _parDrawLevel; // 0(batch), 1(event by event), 2(time-window by time-window), 3(iteration by iteration)
    bool _parDoSmallCurlFinding;
    bool _parReDoFindingNextTimeBin;

    // geometory
    std::vector<double> _parSensorOriginX;
    std::vector<double> _parSensorOriginZ;
    double _parSensorSize;

    // time window
    double _timeWindowMin;
    double _timeWindowMax;
    double _parTimeWindowStep;
    double _parTimeWindowWidth;
    TH1I*  _histRecTime;
    
    // MC information
    std::map<int, const MCParticle*> _mapMCPositron; // key is muon-ID
    std::map<int, const MCParticle*> _mapMCMuon;     // key is muon-ID
    std::map<int,std::vector<int> >  _mapMCTrack;    // <key,value> = <muon-ID, vector of index No>
    

    // Hit information
    std::vector<const RecoHit*> _hitInfoRecoHits;
    std::vector<double>         _hitInfoPhi;
    std::vector<int>            _hitInfoVaneID;
    std::vector<int>            _hitInfoMuonID;
    std::vector<int>            _hitInfoLineNoPhiZ; // indicate which hough-line is clonsed (initial value:-999)
    std::vector<int>            _hitInfoClusterNo;  // cluster number (initial value:-999)
    std::vector<int>            _hitInfoSequenceNo; // sequence number in the cluster (initial value : -999, start from zero)

    // sort
    std::vector<int> _orderVaneID;
    std::vector<int> _orderRecTime;
    std::vector<int> _orderMCTruthTime;
    std::vector<int> _orderZ;
    std::unordered_multimap<int,int> _multimapVaneID;

    // Hough Transformation in HoughTransformPhiZ()
    std::vector<TH2D*>                _houghPhiZHist;
    std::vector<std::vector<double> > _houghPhiZRho;   // rho-value for each hits
    std::vector<double>               _houghPhiZTheta; // step in theta axis

  // Hough Fit in HoughFit_phiz()
    int                    _houghPhiZLineNoMax;
    std::map<int,TF1*>     _houghPhiZFunc; // straight line on phi-z plane. Deactivated line is removed in clearTime();

    // in CalcHoughResidual_phiz()
    std::vector<TH1D*> _houghPhiZHistResi;

    // Clustering
    int _clusterNoMax; // maximum cluster number (initial value:-1)
    std::map<int,std::vector<int> > _mapClusteringResult;      // <key,value> = <clusterNo, vector of index No>, which  is used after reconstructed clusters
    std::map<int,int>               _mapActiveClusterNoLineNo; // <key,value> = <clusterNo, lineNo>, lineNo can be searched from clusterNo.
    std::map<int,int>               _mapActiveClusterNo;       // <key,value> = <clusterNo, end time>, used in Clustering_NextTimeBin
    std::multiset<int>              _setActiveLineNo;          // key is lineNo
    std::set<int>                   _setActiveMuonID;          // key is muon-ID

    // Extrapolation
    int calcPerpLine( double x1,  double y1,  double x2, double y2, double& slope, double& offset ); // calculate perpendicular bisector from the given two points.
    int circleBy3Point ( double x1,  double y1,  double x2, double y2, double x3, double y3,
			 double& x0, double& y0, double& r ); // calculate a circle from the given three points
    int intersectionCircleLine( double x0, double y0, double r, double slope, double offset, double& x1, double& y1, double& x2, double& y2 );
    int extrapolation( int index1, int index2, int index3, int VaneID, double& extrap_r, double& extrap_z,
		       double& x0, double& y0, double& r, double& dphi,
		       bool isneardphi=true );
    // index1 => index2 => index3 => extrapolated point.
    // clockwise and anti-clockwise direction is automatically judged based on total path length from index1 to index3 through index2.
    // the intersection with smaller dphi is returned between two intersections in the case of isneardphi=true.
    // 1(normal), 2(out of r-region), -2(out of z-region), -1(error)
    
    void writeRecoTrks();
    std::vector<const Track*> _tracks;

    void drawEvtDisplay(int nEvent);
    TArc*   _detectorMuonOrbit;
    TLine** _detectorVane;
    TArc*   _detectorPole;
    std::map<int,TMarker*> _mapDecayPointXY;       // key is muon-ID
    std::map<int,TMarker*> _mapDecayPointPhiZ;     // key is muon-ID
    std::map<int,TArrow*>  _mapDecayDirectionXY;   // key is muon-ID
    std::map<int,TArrow*>  _mapDecayDirectionPhiZ; // key is muon-ID
    std::map<int,TArc*>    _mapIdealTrajectory;    // key is muon-ID

    TGraph* _graphHitPointXYAll;
    TGraph* _graphHitPointXYSig;
    TGraph* _graphHitPointXYBkg;
    TGraph* _graphHitPointPhiZAll;
    TGraph* _graphHitPointPhiZSig;
    TGraph* _graphHitPointPhiZBkg;
    double  readMCInfo();
    
    void drawClusteringResult();
    
    std::map<int,TGraph*> _mapClusteredHitXY;   // key is clusterNo
    std::map<int,TGraph*> _mapClusteredHitPhiZ; // kei is clsuterNo

    // for performance study
    bool _indepStudy;
    TTree* _tree;
    std::vector<int>    _performance_nhits;
    std::vector<int>    _performance_ngoodhits;
    std::vector<int>    _performance_ncontinuous_hits;
    std::vector<int>    _performance_muonid;
    std::vector<int>    _performance_clusterno;
    std::vector<int>    _performance_existmuonid;
    std::vector<double> _performance_existenergy;
    std::vector<double> _performance_existpt;
    std::vector<double> _performance_existpz;
    
    int    getMuonIDfromRecoHit     ( const RecoHit* recohit );
    double getMCTruthTimefromRecoHit( const RecoHit* recohit );
    void   checkRecoHits            ( const RecoHit* recohit );
    
    // Time windows
    int   getTimeWindowMin(){ return _timeWindowMin; }
    int   getTimeWindowMax(){ return _timeWindowMax; }
    void  setTimeWindowMin( double time_window_min ){ _timeWindowMin = time_window_min; }
    void  setTimeWindowMax( double time_window_max ){ _timeWindowMax = time_window_max; }
    void  forwardTimeWindow( int nstep=1 ){ _timeWindowMin += nstep*_parTimeWindowStep; _timeWindowMax += nstep*_parTimeWindowStep;  };
    void  makeRecTimeHist( double time_window_start, double time_window_end );
    bool  skipTimeWindow();
    
    void  calcActiveMuonID( double time_window_min=-999, double time_window_max=-999 );
    
    void    makeGraph              ( TGraph* g_all, TGraph* g_sig, TGraph* g_bkg, int sel, double time_window_min=-999, double time_window_max=-999 );
    void    makeGraphXY            ( TGraph* g_all, TGraph* g_sig, TGraph* g_bkg,          double time_window_min=-999, double time_window_max=-999 );
    void    makeGraphPhiZ          ( TGraph* g_all, TGraph* g_sig, TGraph* g_bkg,          double time_window_min=-999, double time_window_max=-999 );

    TGraph* makeGraphClusterNo    ( int sel, int clusterNo=-999, double time_window_min=-999, double time_window_max=-999 );
    TGraph* makeGraphClusterNoXY  (          int clusterNo=-999, double time_window_min=-999, double time_window_max=-999 );
    TGraph* makeGraphClusterNoPhiZ(          int clusterNo=-999, double time_window_min=-999, double time_window_max=-999 );
    TGraph* makeGraphLineNo       ( int sel, int lineNo,         double time_window_min=-999, double time_window_max=-999 );
    TGraph* makeGraphLineNoXY     (          int lineNo,         double time_window_min=-999, double time_window_max=-999 );
    TGraph* makeGraphLineNoPhiZ   (          int lineNo,         double time_window_min=-999, double time_window_max=-999 );
    TGraph* makeGraphMuonID       ( int sel, int muonID,         double time_window_min=-999, double time_window_max=-999 );
    TGraph* makeGraphMuonIDXY     (          int muonID,         double time_window_min=-999, double time_window_max=-999 );
    TGraph* makeGraphMuonIDPhiZ   (          int muonID,         double time_window_min=-999, double time_window_max=-999 );
    void    makeMuonDecayPlot     ();
    int     incrementVaneID( const int vaneID, const int dev=1 );
    int     decrementVaneID( const int vaneID, const int dev=1 );
    double  calcDeltaPhi   ( double phi_start, double phi_end, bool isclockwise );
    bool    isMoveAwayFromOrigin( int fl_incre, int target_vaneID, int origin_vaneID );

    int     getDominantMuonID( std::multiset<int> & cnt_muonid );
  };
}

#endif
