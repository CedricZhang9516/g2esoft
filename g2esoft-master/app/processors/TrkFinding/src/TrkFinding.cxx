// TrkFinding.cxx

#include "TreeProc/TreeManager.h"
#include "TreeProc/LogStream.h"
#include "TrkFinding.h"

#include <iostream>
#include <iomanip>

using namespace g2esoft;
using namespace TreeProc;
using namespace std;

static Factory<Processor, TrkFinding> _f("TrkFinding");

TrkFinding::TrkFinding(const char *name, const char *instName)
  : Processor(name, instName),
    _geomRPole(72),
    _indepStudy(false),
    _histRecTime(nullptr)
{
}

TrkFinding::~TrkFinding(){}

void TrkFinding::config(Parameter *param){
  log("info") << "*****[START: config]*****" << std::endl;

  _parRecoHitName      = param->get("RecoHitName",     string("RecoHits"     ));
  _parTrackName        = param->get("TrackName",       string("Tracks"       ));
  _parMCParticleName   = param->get("MCParticleName",  string("MCParticles"  ));

  _parSensorOriginX = param->get("SensorOriginX",vector<double>{90.96,190.23});
  _parSensorOriginZ = param->get("SensorOriginZ",vector<double>{0.5,99.77});
  _parSensorSize    = param->get("SensorSize",(double)98.77);

  _parGeomNVane     = param->get("NVane",(int)40);

  _parHoughNStepTheta = param->get("HoughNStepTheta", (int) 180);
  _parHoughNStepRho   = param->get("HoughNStepRho",   (int)1000);
  _parHoughRhoMin     = param->get("HoughRhoMin",     (int)-500);
  _parHoughRhoMax     = param->get("HoughRhoMax",     (int) 500);

  _parCutHoughPhiZPeak         = param->get("CutHoughPhiZPeak",         (int) 3);
  _parCutHoughPhiZPeakDeltaRay = param->get("CutHoughPhiZPeakDeltaRay", (int)10);
  _parCutTrackRangeThetaMin    = param->get("CutTrackRangeThetaMin",    (double) 30.0);
  _parCutTrackRangeThetaMax    = param->get("CutTrackRangeThetaMax",    (double)150.0);
  _parCutDeltaRayRangeThetaMax = param->get("CutDeltaRayRangeThetaMax", (double)  3.0);

  _parCutHoughPhiZSeedResi         = param->get("CutHoughPhiZSeedResi",         (double)10.0);
  _parCutHoughPhiZSeedResiDeltaRay = param->get("CutHoughPhiZSeedResiDeltaRay", (double)0.01);

  _parCutExtrapTolerance              = param->get("CutExtrapTolerance",              (double)5.0);
  _parCutExtrapToleranceCoeffDphi     = param->get("CutExtrapToleranceCoeffDphi",     (double)5.0);

  _parCutExtrapMiss    = param->get("CutExtrapMiss",    (int)5);
  _parCutExtrapNoCross = param->get("CutExtrapNoCross", (int)5);

  _parThresholdSuccess = param->get("ThresholdSuccess", (int)5);

  _parTimeWindowStep  = param->get("TimeWindowStep", (double) 5.0);
  _parTimeWindowWidth = param->get("TimeWindowWidth",(double)10.0);

  _parDrawLevel                 = param->get("DrawLevel",(int)0);
  _parDoSmallCurlFinding     = param->get("DoSmallCurlFinding",    (bool)false);
  _parReDoFindingNextTimeBin = param->get("ReDoFindingNextTimeBin",(bool)false);

  _parEvtDispThresholdMomentum = param->get("EvtDispThresholdMomentum",(double)200.0);

  _geomRInner = _parSensorOriginX.at(0);
  _geomROuter = _parSensorOriginX.at(1) + _parSensorSize;
  _geomZMin   = -1.0*(_parSensorOriginZ.at(1)+_parSensorSize);
  _geomZMax   =       _parSensorOriginZ.at(1)+_parSensorSize;

  printParameters();
        
  // register Tracks
  TreeManager::instance()->registerBranch(_parTrackName.c_str(), &_tracks);

  // register Objects for Trking performance study
  if( _indepStudy ){
    _tree = new TTree("evaltrking","evaltrking");
    _tree->Branch("nhits",            &_performance_nhits);
    _tree->Branch("ngoodhits",        &_performance_ngoodhits);
    _tree->Branch("ncontinuous_hits", &_performance_ncontinuous_hits);
    _tree->Branch("muonid",           &_performance_muonid);
    _tree->Branch("clusterno",        &_performance_clusterno);
    _tree->Branch("existmuonid",      &_performance_existmuonid);
    _tree->Branch("existenergy",      &_performance_existenergy);
    _tree->Branch("existpt",          &_performance_existpt);
    _tree->Branch("existpz",          &_performance_existpz);
  }
  if( !_parDrawLevel )gROOT->SetBatch(true);

  _geomPhi = new double[_parGeomNVane];
  for( int ivane=0; ivane<_parGeomNVane; ivane++ ){
    _geomPhi[ivane] = -2.0*TMath::Pi()*ivane/(double)_parGeomNVane-1.e-9;
    log("debug") << "ivane = " << ivane << ", phi = " << _geomPhi[ivane] << std::endl;
  }

  _canvas1 = new TCanvas( "canvas1", "canvas1", 1200, 800 );
  _canvas1->Divide(3,2);
  _canvas1->Draw();

  _detectorMuonOrbit = new TArc( 0, 0, 330 ); // hard-code
  _detectorMuonOrbit->SetFillStyle(0);
  _detectorMuonOrbit->SetLineColor(kGray);

  // Detector Geometory
  _detectorVane = new TLine*[_parGeomNVane];
  for( unsigned int ivane=0; ivane < _parGeomNVane; ivane++ ){
    double phi = _geomPhi[ivane];
    _detectorVane[ivane] = new TLine( _geomRInner*TMath::Cos(phi),
				      _geomRInner*TMath::Sin(phi),
				      _geomROuter*TMath::Cos(phi),
				      _geomROuter*TMath::Sin(phi)
				      );
    _detectorVane[ivane]->SetLineWidth(1);
    _detectorVane[ivane]->SetLineColor(kGray);
  }
	 
  _detectorPole = new TArc( 0, 0, _geomRPole );
  _detectorPole->SetFillStyle(0);
  _detectorPole->SetLineColor(kGray);

  _graphHitPointXYAll   = new TGraph();
  _graphHitPointXYSig   = new TGraph();
  _graphHitPointXYBkg   = new TGraph();
  _graphHitPointXYAll  ->SetMarkerStyle(20);
  _graphHitPointXYSig  ->SetMarkerStyle(20);
  _graphHitPointXYBkg  ->SetMarkerStyle(20);
  _graphHitPointXYAll  ->SetMarkerSize(0.5);
  _graphHitPointXYSig  ->SetMarkerSize(0.5);
  _graphHitPointXYBkg  ->SetMarkerSize(0.5);
  _graphHitPointXYAll  ->SetMarkerColor(1);
  _graphHitPointXYSig  ->SetMarkerColor(2);
  _graphHitPointXYBkg  ->SetMarkerColor(1);

  _graphHitPointPhiZAll = new TGraph();
  _graphHitPointPhiZSig = new TGraph();
  _graphHitPointPhiZBkg = new TGraph();
  _graphHitPointPhiZAll->SetMarkerStyle(20);
  _graphHitPointPhiZSig->SetMarkerStyle(20);
  _graphHitPointPhiZBkg->SetMarkerStyle(20);
  _graphHitPointPhiZAll->SetMarkerSize(0.5);
  _graphHitPointPhiZSig->SetMarkerSize(0.5);
  _graphHitPointPhiZBkg->SetMarkerSize(0.5);
  _graphHitPointPhiZAll->SetMarkerColor(1);
  _graphHitPointPhiZSig->SetMarkerColor(2);
  _graphHitPointPhiZBkg->SetMarkerColor(1);

  return;
}

void TrkFinding::process(int nEvent, Parameter* ){
  log("info") << "start TrkFinding process for event: " << nEvent << endl;

  // clear event;
  clearEvent();

  // input: recoHits
  const vector<const RecoHit*> &recoHits = TreeManager::instance()->getBranchVec<RecoHit>(_parRecoHitName.c_str());
  if( recoHits.size()==0 ) return;

  // input variables
  log("debug") << "input variables" << std::endl;
  for( unsigned int ivec=0; ivec<recoHits.size(); ivec++ ) inputHits( ivec, recoHits.at(ivec) );
  
  // calculation order
  calcOrder();

  // print hit information
  //printHitInfoVaneIDOrder(); // printHitInfo(), printHitInfoVaneIDOrder(), printHitInfoRecTimeOrder(), printHitInfoMCTruthTimeOrder(), printHitInfoZOrder();

  // start clustering process
  readMCInfo();
  makeMuonDecayPlot();

  log("debug") << "start track finding algorithm !" << std::endl;

  int     cnt_cycle         =  0;
  int     fl_cycle_end      =  0;
  double  time_window_start = ( _hitInfoRecoHits.size() ? _hitInfoRecoHits.at(_orderRecTime.at(0                 ))->_time : -999 );
  double  time_window_end   = ( _hitInfoRecoHits.size() ? _hitInfoRecoHits.at(_orderRecTime.at(_hitInfoRecoHits.size()-1))->_time : -999 );

  makeRecTimeHist(time_window_start, time_window_end);

  bool fl_next_time      = false;
  bool fl_first_time     = true;
  
  setTimeWindowMin( time_window_start );
  setTimeWindowMax( time_window_start + _parTimeWindowWidth );
  calcActiveMuonID();
  log("info") << "[Time-Window] " << _timeWindowMin << " ~ " << _timeWindowMax << std::endl;
  int cnt_show_time = 0;


  while(1){ // START-WHILE
    if( fl_first_time ){
	;
    }else if( fl_next_time ){ // move next time bin
      clearTime();
      if( _timeWindowMax > time_window_end ) break;
      
      forwardTimeWindow();
      if( skipTimeWindow() ) continue;
      log("info") << "[Time-Window] " << _timeWindowMin << " ~ " << _timeWindowMax << std::endl;
      calcActiveMuonID();

      // Extrapolation of hits at next time bin from already reconstructed clusters
      if( _parReDoFindingNextTimeBin ) clusteringNextTimeBin(_timeWindowMax-_parTimeWindowStep, _timeWindowMax); // to be checked

      cnt_cycle = 0;
      cnt_show_time++;

    }

    if( fl_next_time || fl_first_time ) drawEvtDisplay(nEvent);
    if( fl_next_time  ) fl_next_time  = false;
    if( fl_first_time ) fl_first_time = false;

    fl_cycle_end = 0;

    // Hough Transformation (phi-Z)
    houghTransformPhiZ();

    // Search Hough-fit line
    if( houghFitPhiZ() < 0 ) fl_cycle_end = 1;

    // Find hit points close to Hough-fit line
    if( !fl_cycle_end && calcHoughResidualPhiZ() < 0 ) fl_cycle_end = 1;

    //clustering
    if( !fl_cycle_end && clusteringShort() < 0 ) fl_cycle_end = 1;

    drawClusteringResult();
    if( _parDrawLevel>2 ){
      _canvas1->Update();
      _canvas1->WaitPrimitive();
    }

    cnt_cycle++;
    log("debug") << "One cycle of track finding has been finished : Cycle#" << cnt_cycle << std::endl << std::endl;
    if( fl_cycle_end ){
      fl_next_time = true;
      log("debug") << "Track finding for one time window has been completely finished and we will move to next time bin : Cycle#" << cnt_cycle << std::endl << std::endl;
      if( _parDrawLevel>1 ){
	_canvas1->Update();
	//printHitInfoVaneIDOrder();
	_canvas1->WaitPrimitive();
      }
    }

  } // END-WHILE

  if( _parDrawLevel ) printHitInfoVaneIDOrder();
  // Evaluate track finding performance
  //std::vector<std::pair<int,int> >* result = new std::vector<pair<int,int> >;
  //judgeFinding(result);
  evaluateTrkFinding();

  // Write output
  writeRecoTrks();

  if( _parDrawLevel>0 ){
    _canvas1->Update();
    _canvas1->WaitPrimitive();
  }

  return;  
}


void TrkFinding::finish(Parameter *){
  log("debug") << "*****[START: finish]*****" << std::endl;
  clearEvent();

  delete _detectorMuonOrbit;
  for( int ivane=0; ivane<_parGeomNVane; ivane++ ) delete _detectorVane[ivane];
  delete[] _detectorVane;
  delete _detectorPole;

  delete _graphHitPointXYAll;
  delete _graphHitPointXYSig;
  delete _graphHitPointXYBkg;
  delete _graphHitPointPhiZAll;
  delete _graphHitPointPhiZSig;
  delete _graphHitPointPhiZBkg;

  std::map<int,TMarker*>::iterator it_decpoint_xy = _mapDecayPointXY.begin();
  while( it_decpoint_xy != _mapDecayPointXY.end() ){
    delete it_decpoint_xy->second;
    _mapDecayPointXY.erase(it_decpoint_xy++);
  }

  std::map<int,TMarker*>::iterator it_decpoint_phiz = _mapDecayPointPhiZ.begin();
  while( it_decpoint_phiz != _mapDecayPointPhiZ.end() ){
    delete it_decpoint_phiz->second;
    _mapDecayPointPhiZ.erase(it_decpoint_phiz++);
  }

  std::map<int,TArrow*>::iterator it_decvec_xy = _mapDecayDirectionXY.begin();
  while( it_decvec_xy != _mapDecayDirectionXY.end() ){
    delete it_decvec_xy->second;
    _mapDecayDirectionXY.erase(it_decvec_xy++);
  }

  std::map<int,TArrow*>::iterator it_decvec_phiz = _mapDecayDirectionPhiZ.begin();
  while( it_decvec_phiz != _mapDecayDirectionPhiZ.end() ){
    delete it_decvec_phiz->second;
    _mapDecayDirectionPhiZ.erase(it_decvec_phiz++);
  }

  std::map<int,TArc*>::iterator it_ideal_trk = _mapIdealTrajectory.begin();
  while( it_ideal_trk != _mapIdealTrajectory.end() ){
    delete it_ideal_trk->second;
    _mapIdealTrajectory.erase(it_ideal_trk++);
  }
  
  if( _indepStudy ) _tree->Write();
  
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TrkFinding::inputHits( int index, const RecoHit* recohit ){
  _hitInfoRecoHits.push_back         ( recohit                                      );
  _hitInfoPhi.push_back              ( calcPhi(recohit->_pos.Y(),recohit->_pos.X()) );
  _hitInfoVaneID.push_back           ( getVaneID(_hitInfoPhi.at(_hitInfoPhi.size()-1))            );
  _hitInfoMuonID.push_back           ( getMuonIDfromRecoHit(recohit)                );
  _hitInfoLineNoPhiZ.push_back ( -999                                         );
  _hitInfoClusterNo.push_back        ( -999                                         );
  _hitInfoSequenceNo.push_back       ( -999                                         );

  return;
}

void TrkFinding::clearTime( ){
  log("debug") << "*****[START: clearTime]*****" << std::endl;
  // HoughTransform_phiz()
  for( unsigned int ivec=0; ivec<_houghPhiZHist.size(); ivec++ ) delete _houghPhiZHist.at(ivec);
  _houghPhiZHist.clear();
  _houghPhiZRho.clear();
  _houghPhiZTheta.clear();

  // HoughFit_phiz()

  // deactivate clusterNo
  std::map<int,int>::iterator it_cluster = _mapActiveClusterNo.begin();
  Int_t end_time, clusterNo, lineNo;
  log("debug") << "Number of active clusters : " << getNActiveClusters() << std::endl;
  while( it_cluster != _mapActiveClusterNo.end() ){
    end_time  = it_cluster->second;
    clusterNo = it_cluster->first;
    lineNo    = _mapActiveClusterNoLineNo.at(clusterNo);
    log("debug") << "    clusterNo = "     << clusterNo
		<< ", lineNo = "          << lineNo
		<< ", end time = "        << end_time
		<< ", time_window_min = " << _timeWindowMin
		<< ", time_window_max = " << _timeWindowMax
		<< std::endl;
    if( end_time < _timeWindowMax - _parTimeWindowStep ){ // condition for deactivation (it may be better to adjust)
      _mapActiveClusterNoLineNo.erase(clusterNo);    // deactivate lineNo
      _mapActiveClusterNo.erase       (it_cluster++); // deactivate clusterNo
      _setActiveLineNo.erase          (_setActiveLineNo.find(lineNo));
      _mapClusteredHitXY.erase       (clusterNo);
      _mapClusteredHitPhiZ.erase     (clusterNo);
      log("debug") << "  ===> deactivate clusterNo = " << clusterNo << ", lineNo = " << lineNo << std::endl;
    }else ++it_cluster;
  }

  // deactivate lineNo
  std::map<int,TF1*>::iterator it_houghPhiZFunc = _houghPhiZFunc.begin();
  log("debug") << "Number of active lines : " << getNLines() << std::endl;
  while( it_houghPhiZFunc != _houghPhiZFunc.end() ){
    int  lineNo = it_houghPhiZFunc->first;
    TF1* func   = it_houghPhiZFunc->second;
    if( _setActiveLineNo.count(lineNo)==0 ){      
      delete func;
      _houghPhiZFunc.erase(it_houghPhiZFunc++);
      log("debug") << "    ===> deactivate lineNo = " << lineNo << std::endl;
    }else{
      log("debug") << "    lineNo = " << lineNo << " (active)" << std::endl;
      ++it_houghPhiZFunc;
    }
  }
  
  log("debug") << "Nclsuter(active) = " << getNActiveClusters() << ", "
	      << "Nlines  (active) = " << getNActiveLines()    << " (after deactivation)" << std::endl;
  
  // CalcHoughResidual_phiz()
  for( unsigned ivec=0; ivec<_houghPhiZHistResi.size(); ivec++ ) delete _houghPhiZHistResi.at(ivec);
  _houghPhiZHistResi.clear();

  return;
}


void TrkFinding::clearEvent(){
  log("debug") << "*****[START : clearEvent]*****" << std::endl;
  clearTime();

  for( unsigned int itrk=0; itrk<_tracks.size(); itrk++ ) delete _tracks.at(itrk);
  _tracks.clear();

  // MC information
  _mapMCMuon.clear();
  _mapMCPositron.clear();
  std::map<int, std::vector<int> >::iterator it_MCTrack = _mapMCTrack.begin();
  while( it_MCTrack != _mapMCTrack.end() ){
    it_MCTrack->second.clear();
    ++it_MCTrack;
  }
  _mapMCTrack.clear();  
  

  // Hit information
  _hitInfoRecoHits.clear();
  _hitInfoPhi.clear();
  _hitInfoVaneID.clear();
  _hitInfoMuonID.clear();
  _hitInfoLineNoPhiZ.clear();
  _hitInfoClusterNo.clear();
  _hitInfoSequenceNo.clear();

  // order
  _orderVaneID.clear();
  _orderRecTime.clear();
  _orderMCTruthTime.clear();
  _orderZ.clear();
  _multimapVaneID.clear();

  // HoughTransform_phiz()

  // HoughFit_phiz()
  std::map<int,TF1*>::iterator it_houghPhiZFunc = _houghPhiZFunc.begin();
  while( it_houghPhiZFunc != _houghPhiZFunc.end() ){
    delete it_houghPhiZFunc->second;
    _houghPhiZFunc.erase(it_houghPhiZFunc++);
  }
  _houghPhiZFunc.clear();
  
  _houghPhiZLineNoMax = -1;

  // clustering
  _clusterNoMax = -1;
  _mapClusteringResult.clear();
  _mapActiveClusterNo.clear();
  _mapActiveClusterNoLineNo.clear();
  _setActiveLineNo.clear();

  std::map<int,TGraph*>::iterator it_mapClusteredHitXY = _mapClusteredHitXY.begin();
  while( it_mapClusteredHitXY != _mapClusteredHitXY.end() ){
    delete it_mapClusteredHitXY->second;
    _mapClusteredHitXY.erase(it_mapClusteredHitXY++);
  }
  
  std::map<int,TGraph*>::iterator it_mapClusteredHitPhiZ = _mapClusteredHitPhiZ.begin();
  while( it_mapClusteredHitPhiZ != _mapClusteredHitPhiZ.end() ){
    delete it_mapClusteredHitPhiZ->second;
    _mapClusteredHitPhiZ.erase(it_mapClusteredHitPhiZ++);
  }
  
  delete _histRecTime;
  _histRecTime = nullptr;

  return;
}

unsigned int TrkFinding::getVaneID( double phi ){
  while(phi<0.) phi += 2*TMath::Pi();
  while(phi>2*TMath::Pi()) phi -= 2*TMath::Pi();
  int vaneID = (int)((1.-phi/(2.*TMath::Pi()))*_parGeomNVane+0.5);
  vaneID = vaneID%_parGeomNVane;
  
  return vaneID;
}

double TrkFinding::calcPhi( double y, double x ){
  double phi = TMath::ATan2(y,x);
  if( phi<0 ) phi += 2*TMath::Pi();
  return phi;
}

void TrkFinding::printParameters(){
  log("info") << "*****[START : printParameter]**************************"     << std::endl;
  log("info") //<< "SimHitName = "                          << _parSimHitName.c_str()                  << std::endl
	      //<< "StripClusterName = "                    << _parStripClusterName.c_str()            << std::endl
	      << "RecoHitName = "                         << _parRecoHitName.c_str()                 << std::endl
	      << "TrackName = "                           << _parTrackName.c_str()                   << std::endl
	      << "MCParticleName = "                      << _parMCParticleName.c_str()              << std::endl
	      << "NVane  = "                              << _parGeomNVane                           << std::endl
	      << "geomRInner = "                          << _geomRInner                             << std::endl
	      << "geomROuter = "                          << _geomROuter                             << std::endl
	      << "geomZMin = "                            << _geomZMin                               << std::endl
	      << "geomZMax = "                            << _geomZMax                               << std::endl
	      << "HoughNStepTheta= "                      << _parHoughNStepTheta                     << std::endl
	      << "HoughNStepRho = "                       << _parHoughNStepRho                       << std::endl
	      << "HoughRhoMin = "                         << _parHoughRhoMin                         << std::endl
	      << "HoughRhoMax = "                         << _parHoughRhoMax                         << std::endl
	      << "CutHoughPhiZPeak = "                    << _parCutHoughPhiZPeak                    << std::endl
	      << "CutHoughPhiZPeakDeltaRay = "            << _parCutHoughPhiZPeakDeltaRay            << std::endl
	      << "CutTrackRangeThetaMin = "               << _parCutTrackRangeThetaMin               << std::endl
	      << "CutTrackRangeThetaMax = "               << _parCutTrackRangeThetaMax               << std::endl
	      << "CutDeltaRayRangeThetaMax = "            << _parCutDeltaRayRangeThetaMax            << std::endl
	      << "CutHoughPhiZSeedResi = "         << _parCutHoughPhiZSeedResi         << std::endl
	      << "CutHoughPhiZSeedResiDeltaRay = " << _parCutHoughPhiZSeedResiDeltaRay << std::endl
	      << "CutExtrapTolerance = "                  << _parCutExtrapTolerance                  << std::endl
	      << "CutExtrapToleranceCoeffDphi = "         << _parCutExtrapToleranceCoeffDphi         << std::endl
              << "CutExtrapMiss = "                       << _parCutExtrapMiss                       << std::endl
	      << "CutExtrapNoCross = "                    << _parCutExtrapNoCross                    << std::endl
	      << "ThresholdSuccess = "                    << _parThresholdSuccess                    << std::endl
	      << "TimeWindowStep = "                      << _parTimeWindowStep                      << std::endl
	      << "TimeWindowWidth = "                     << _parTimeWindowWidth                     << std::endl
	      << "DrawLevel = "                           << _parDrawLevel                           << std::endl
	      << "DoSmallCurlFinding = "                  << _parDoSmallCurlFinding                  << std::endl
	      << "ReDoFindingNextTimeBin = "              << _parReDoFindingNextTimeBin              << std::endl
    //<< "EvtDispThresholdSignalEnergy = "        << _parEvtDispThresholdSignalEnergy        << std::endl
	      << "EvtDispThresholdMomentum = "            << _parEvtDispThresholdMomentum            << std::endl
    ;
    
  return;
}

void TrkFinding::printHitInfo(int mode){
  log("debug") << "*****[START : printHitInfo]*****" << std::endl;
  if     ( mode==0 ) log("debug") <<                     "Nvec = " << _hitInfoRecoHits.size()          << std::endl;
  else if( mode==1 ) log("debug") <<      "VaneID-Order : Nvec = " << _orderVaneID.size()      << std::endl;
  else if( mode==2 ) log("debug") <<     "RecTime-Order : Nvec = " << _orderRecTime.size()     << std::endl;
  else if( mode==3 ) log("debug") << "MCTruthTime-Order : Nvec = " << _orderMCTruthTime.size() << std::endl;
  else if( mode==4 ) log("debug") <<           "Z-Order : Nvec = " << _orderZ.size()           << std::endl;
  else{
    log("debug") << "[ERROR] Wrong mode in printHitInfo()" << std::endl;
    return;
  }
  int index;
  for( int ivec=0; ivec<_hitInfoRecoHits.size(); ivec++ ){
    if     ( mode==0 ) index = ivec;
    else if( mode==1 ) index = _orderVaneID.at(ivec);
    else if( mode==2 ) index = _orderRecTime.at(ivec);
    else if( mode==3 ) index = _orderMCTruthTime.at(ivec);
    else if( mode==4 ) index = _orderZ.at(ivec);
    log("debug") << "     "
		<< setw(3) << right << ivec << " : "
		<< setw(3) << right << index << " : (x,y,z) = ("
		<< setw(7) << right << Form("%.2f",_hitInfoRecoHits.at(index)->_pos.X()   ) << ", "
		<< setw(7) << right << Form("%.2f",_hitInfoRecoHits.at(index)->_pos.Y()   ) << ", "
		<< setw(7) << right << Form("%.2f",_hitInfoRecoHits.at(index)->_pos.Z()   ) << "), (R,phi) = ("
		<< setw(7) << right << Form("%.2f",_hitInfoRecoHits.at(index)->_pos.Perp()) << ", "
		<< setw(5) << right << Form("%.2f",_hitInfoPhi.at(index)                  ) << "), vane-ID = "
		<< setw(2) << right << Form("%d",  _hitInfoVaneID.at(index)               ) << ", rec-time = "
		<< setw(6) << right << Form("%.2f",_hitInfoRecoHits.at(index)->_time      ) << ", mc-time = "
		<< setw(8) << right << Form("%.3f",getMCTruthTimefromRecoHit(_hitInfoRecoHits.at(ivec))) << ", MuonID = "
		<< setw(3) << right << Form("%d",  _hitInfoMuonID.at(index)          ) << ", Type = "
		<< setw(1) << right << Form("%d",  _hitInfoRecoHits.at(index)->_trueType ) << " : "
		<< setw(4) << right << _hitInfoLineNoPhiZ.at(index)              << ", "
		<< setw(4) << right << _hitInfoClusterNo.at       (index)              << ", "
		<< setw(4) << right << _hitInfoSequenceNo.at      (index)
		<< endl;
  }
  return;
}

void TrkFinding::calcOrder(){
  log("debug") << "*****[START: calcOrder]*****" << std::endl;
  if( _orderVaneID.size()      ) _orderVaneID.clear();
  if( _orderRecTime.size()     ) _orderRecTime.clear();
  if( _orderMCTruthTime.size() ) _orderMCTruthTime.clear();
  if( _orderZ.size()           ) _orderZ.clear();
  if( _multimapVaneID.size()   ) _multimapVaneID.clear();

  multimap<int,   double> tMap_VaneID;
  multimap<double,double> tMap_RecTime;
  multimap<double,double> tMap_MCTruthTime;
  multimap<double,double> tMap_Z;

  for( unsigned int ivec=0; ivec< _hitInfoRecoHits.size(); ivec++ ){
    tMap_VaneID.insert     ( make_pair(_hitInfoVaneID.at(ivec),                             ivec) );
    tMap_RecTime.insert    ( make_pair(_hitInfoRecoHits.at(ivec)->_time,                    ivec) );
    tMap_MCTruthTime.insert( make_pair(getMCTruthTimefromRecoHit(_hitInfoRecoHits.at(ivec)),ivec) );
    tMap_Z.insert          ( make_pair(_hitInfoRecoHits.at(ivec)->_pos.Z(),                 ivec) );
    _multimapVaneID.insert                  ( make_pair(_hitInfoVaneID.at(ivec), ivec) );
  }

  multimap<int,   double>::iterator it_VaneID      = tMap_VaneID.begin();
  multimap<double,double>::iterator it_RecTime     = tMap_RecTime.begin();
  multimap<double,double>::iterator it_MCTruthTime = tMap_MCTruthTime.begin();
  multimap<double,double>::iterator it_Z           = tMap_Z.begin();
  while( it_VaneID      != tMap_VaneID.end()      ){ _orderVaneID.push_back     ( (*it_VaneID     ).second ); it_VaneID++;      }
  while( it_RecTime     != tMap_RecTime.end()     ){ _orderRecTime.push_back    ( (*it_RecTime    ).second ); it_RecTime++;     }
  while( it_MCTruthTime != tMap_MCTruthTime.end() ){ _orderMCTruthTime.push_back( (*it_MCTruthTime).second ); it_MCTruthTime++; }
  while( it_Z           != tMap_Z.end()           ){ _orderZ.push_back          ( (*it_Z          ).second ); it_Z++;           }

  return;
}


int TrkFinding::calcPerpLine( double x1, double y1, double x2, double y2, double& slope, double& offset ){
  if( TMath::Abs(x1-x2)<1.0e-5 && TMath::Abs(y1-y2)<1.0e-5 ){
    //std::cerr << "[ABORT] Same points are input" << std::endl;
    return -1;
  }
  double slope1 = ( x2!=x1 ? (y2-y1)/(x2-x1) : 0.0 );
  if( TMath::Abs(slope1)<1.0e-5 ){
    //std::cerr << "[ABORT] Invalid slope : " << slope1 << std::endl;
    return -1;
  }
  slope = -1.0/slope1;
  offset = 0.5*((y1+y2) - slope*(x1+x2));
  return 1;
}

int TrkFinding::circleBy3Point( double x1, double y1, double x2, double y2, double x3, double y3,
				double& x0, double& y0, double& r ){
  if( (TMath::Abs(x1-x2)<1.0e-5 && TMath::Abs(y1-y2)<1.0e-5) ||
      (TMath::Abs(x2-x3)<1.0e-5 && TMath::Abs(y2-y3)<1.0e-5) ||
      (TMath::Abs(x3-x1)<1.0e-5 && TMath::Abs(y3-y1)<1.0e-5) ){
    //std::cerr << "[ABORT] Same points are input" << std::endl;
    return -1;
  }

  double slope1;
  double slope2;
  double offset1;
  double offset2;

  calcPerpLine( x1, y1, x2, y2, slope1, offset1 );
  calcPerpLine( x2, y2, x3, y3, slope2, offset2 );
  if( TMath::Abs(slope2-slope1)<1.0e-5 ){
    //std::cerr << "[ABORT] Three points on straight line" << std::endl;
    return -1;
  }
  x0 = (offset2-offset1)/(slope1-slope2);
  y0 = slope1 * x0 + offset1;
  r = sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) );
  return 1;
}

int TrkFinding::intersectionCircleLine( double x0, double y0, double r, double slope, double offset, double& x1, double& y1, double& x2, double& y2 ){
  // circle : (x-x0)^2 + (y-y0)^2 = r^2
  // line   : y = slope*x + offset
  // c2*x^2 + c1*x^1 + c0 = 0;

  double c2 = slope*slope + 1.0;
  double c1 = (offset-y0)*slope - x0;
  double c0 = x0*x0 + (offset-y0)*(offset-y0) - r*r;
  if( c1*c1 - c2*c0 < 0 ) return -1;

  // Two intersections: (x1 is always larger than x2)
  double D = sqrt(c1*c1-c2*c0);
  x1 = (-c1 + D)/c2;
  x2 = (-c1 - D)/c2;
  y1 = slope*x1 + offset;
  y2 = slope*x2 + offset;
  
  return 1;

}


int TrkFinding::extrapolation( int index1, int index2, int index3, int vaneID, double& extrap_r, double& extrap_z, double& x0, double& y0, double& r, double& dphi, bool isneardphi ){
  log("debug") << "*****[START: extrapolation]*****" << std::endl;
  if( circleBy3Point( _hitInfoRecoHits.at(index1)->_pos.X(),
		      _hitInfoRecoHits.at(index1)->_pos.Y(),
		      _hitInfoRecoHits.at(index2)->_pos.X(),
		      _hitInfoRecoHits.at(index2)->_pos.Y(),
		      _hitInfoRecoHits.at(index3)->_pos.X(),
		      _hitInfoRecoHits.at(index3)->_pos.Y(),
		      x0, y0, r )<0 ){
    dphi = 999;
    return -1; // derive circle(x0,y0,r) from given three points.
  }

  double slope_vane = TMath::Tan(_geomPhi[vaneID]);
  
  double phi1 = calcPhi( _hitInfoRecoHits.at(index1)->_pos.Y()-y0, _hitInfoRecoHits.at(index1)->_pos.X()-x0 );
  double phi2 = calcPhi( _hitInfoRecoHits.at(index2)->_pos.Y()-y0, _hitInfoRecoHits.at(index2)->_pos.X()-x0 );
  double phi3 = calcPhi( _hitInfoRecoHits.at(index3)->_pos.Y()-y0, _hitInfoRecoHits.at(index3)->_pos.X()-x0 );
  double dphi_clockwise = TMath::Abs(
				     calcDeltaPhi(phi1,phi2,true)+
				     calcDeltaPhi(phi2,phi3,true)
				     );
  double dphi_anticlockwise = TMath::Abs(
					 calcDeltaPhi(phi1,phi2,false)+
					 calcDeltaPhi(phi2,phi3,false)
					   );
  bool fl_clockwise; // clockwise and anti-clockwise direction is automatically judged based on total path length from index1 to index3 through index2.
  if( dphi_clockwise < dphi_anticlockwise ) fl_clockwise = true;  
  else                                      fl_clockwise = false;

  double x_p,y_p;
  double x_n,y_n;
  int fl_answer = intersectionCircleLine( x0, y0, r, slope_vane, 0.0, x_p, y_p, x_n, y_n );
  if( fl_answer<0 ){
    dphi = 999;
    return -1;
  }
  double phi_p = calcPhi( y_p-y0, x_p-x0 );
  double phi_n = calcPhi( y_n-y0, x_n-x0 );

  log("debug") << "     Extrapolated Circle : x0 = " << x0 << ", y0 = " << y0 << ", r  = " << r << " : fl_clockwise = " << fl_clockwise << std::endl
	       <<   "x1 = "           << _hitInfoRecoHits.at(index1)->_pos.X()
	       << ", y1 = "           << _hitInfoRecoHits.at(index1)->_pos.Y()
	       << ", z1 = "           << _hitInfoRecoHits.at(index1)->_pos.Z()
	       << ", phi1 = "         << phi1
	       << " : Vane-ID = "     << _hitInfoVaneID.at(index1) << std::endl
	       <<   "x2 = "           << _hitInfoRecoHits.at(index2)->_pos.X()
	       << ", y2 = "           << _hitInfoRecoHits.at(index2)->_pos.Y()
	       << ", z2 = "           << _hitInfoRecoHits.at(index2)->_pos.Z()
	       << ", phi2 = "         << phi2
	       << " : Vane-ID = "     << _hitInfoVaneID.at(index2) << std::endl
	       <<   "x3 = "           << _hitInfoRecoHits.at(index3)->_pos.X()
	       << ", y3 = "           << _hitInfoRecoHits.at(index3)->_pos.Y()
	       << ", z3 = "           << _hitInfoRecoHits.at(index3)->_pos.Z()
	       << ", phi3 = "         << phi3
	       << " : Vane-ID = "     << _hitInfoVaneID.at(index3) << std::endl
	       << "Target Vane-ID = " << vaneID << std::endl
	       <<   "x_p = "          << x_p
	       << ", y_p = "          << y_p
	       << ", phi_p = "        << phi_p << std::endl
	       <<   "x_n = "          << x_n
	       << ", y_n = "          << y_n
	       << ", phi_n = "        << phi_n << std::endl;

  double extrap_x;
  double extrap_y;
  double extrap_phi;

  double dphi_p = calcDeltaPhi(phi3, phi_p,fl_clockwise);
  double dphi_n = calcDeltaPhi(phi3, phi_n,fl_clockwise);

  if( ( isneardphi && TMath::Abs(dphi_p) >  TMath::Abs(dphi_n)) ||
      (!isneardphi && TMath::Abs(dphi_p) <= TMath::Abs(dphi_n)) ){
    extrap_x   = x_n;
    extrap_y   = y_n;
    extrap_phi = phi_n;
  }else{
    extrap_x   = x_p;
    extrap_y   = y_p;
    extrap_phi = phi_p;
  }
  
  if( (cos(_geomPhi[vaneID])> 0. && extrap_x <=0.0) ||
      (cos(_geomPhi[vaneID])<=0. && extrap_x >=0.0)
      ){
    dphi = 999;
    return -1;
  }

  dphi = calcDeltaPhi( phi3, extrap_phi,fl_clockwise);

  extrap_z = _hitInfoRecoHits.at(index3)->_pos.Z() + (_hitInfoRecoHits.at(index3)->_pos.Z()-_hitInfoRecoHits.at(index1)->_pos.Z())/calcDeltaPhi(phi1,phi3,fl_clockwise)*calcDeltaPhi(phi3,extrap_phi,fl_clockwise);  
  extrap_r = sqrt(extrap_x*extrap_x + extrap_y*extrap_y);
  log("debug") << "     extrap_x = " << extrap_x << ", extrap_y = " << extrap_y << ", extrp_r = " << extrap_r << ", extrap_z = " << extrap_z << ", extrap_phi = " << extrap_phi << ", dphi = " << dphi << std::endl;

  if( extrap_z < _geomZMin   || extrap_z > _geomZMax   ){
    dphi = 999;
    return -2;
  }
  if( extrap_r < _geomRInner || extrap_r > _geomROuter ) return  2; 
  
  return 1;
}

void TrkFinding::houghTransformPhiZ(){
  log("debug") << "*****[START: houghTransformPhiZ]*****" << std::endl;
  TH2D** hist_hough_phiz_vane = new TH2D*[_parGeomNVane];

  for( unsigned int ivane=0; ivane<_parGeomNVane; ivane++ ){
    hist_hough_phiz_vane[ivane] = new TH2D( Form("hist_hough_phiz_vane_%d_%d",                                       _houghPhiZHist.size(),ivane),
					    Form("Hough%d_%d(#phi-Z#rightarrow#theta-#rho);#theta [#circ];#rho [mm]",_houghPhiZHist.size(),ivane),
					    _parHoughNStepTheta, 0, 180, _parHoughNStepRho, _parHoughRhoMin, _parHoughRhoMax );
  }
  TH2D* hist_hough_phiz = new TH2D( Form("hist_hough_phiz_%d",                                            _houghPhiZHist.size()),
				    Form("Hough%d(#phi-Z#rightarrow#theta-#rho);#theta [#circ];#rho [mm]",_houghPhiZHist.size()),
				    _parHoughNStepTheta, 0, 180, _parHoughNStepRho, _parHoughRhoMin, _parHoughRhoMax );
  _houghPhiZRho.clear();
  _houghPhiZTheta.clear();

  const double deg2rad  = TMath::Pi()/180.;
  const double rad2deg  = 180./TMath::Pi();
  const double step2deg = 180./(double)_parHoughNStepTheta;
  double theta, rho;

  for( unsigned int ivec=0; ivec<_hitInfoRecoHits.size(); ivec++ ){
    std::vector<double> tmp_vector;
    if( _timeWindowMin > _hitInfoRecoHits.at(ivec)->_time || _timeWindowMax < _hitInfoRecoHits.at(ivec)->_time ) continue; // time-window is "min <= time <= max"
    if( _hitInfoClusterNo.at(ivec) >= -1 ) continue; // to be checked : -1 or 0
    
    for( unsigned int istep=0; istep<_parHoughNStepTheta; istep++ ){
      theta = istep*step2deg;
      rho = _hitInfoPhi.at(ivec)*rad2deg*TMath::Cos(theta*deg2rad)+_hitInfoRecoHits.at(ivec)->_pos.Z()*TMath::Sin(theta*deg2rad);
      tmp_vector.push_back( rho );
      if( (_parCutDeltaRayRangeThetaMax < theta && theta < _parCutTrackRangeThetaMin) || (_parCutTrackRangeThetaMax < theta) ) continue;
      hist_hough_phiz_vane[_hitInfoVaneID.at(ivec)]->SetBinContent( hist_hough_phiz_vane[_hitInfoVaneID.at(ivec)]->FindBin(theta+1.0e-5, rho), 1.0 );
    }
    _houghPhiZRho.push_back( tmp_vector );
  }
  for( unsigned int istep=0; istep<_parHoughNStepTheta; istep++ ) _houghPhiZTheta.push_back(istep/180.0*_parHoughNStepTheta);

  for( unsigned int ivane=0; ivane<_parGeomNVane; ivane++ ) hist_hough_phiz->Add( hist_hough_phiz_vane[ivane] );
  _houghPhiZHist.push_back( hist_hough_phiz );
  for( unsigned int ivane=0; ivane<_parGeomNVane; ivane++ ) delete hist_hough_phiz_vane[ivane];

  return;
}


int TrkFinding::houghFitPhiZ(){
  log("debug") << "*****[START : houghFitPhiZ]*****" << std::endl;
  // Search Line
  int max_xbin, max_ybin, max_zbin;
  // Search for standard track first
  _houghPhiZHist.at(_houghPhiZHist.size()-1)->GetXaxis()->SetRangeUser(_parCutTrackRangeThetaMin, _parCutTrackRangeThetaMax);
  _houghPhiZHist.at(_houghPhiZHist.size()-1)->GetMaximumBin(max_xbin, max_ybin, max_zbin);
  double theta = _houghPhiZHist.at(_houghPhiZHist.size()-1)->GetXaxis()->GetBinLowEdge(max_xbin);
  double rho   = _houghPhiZHist.at(_houghPhiZHist.size()-1)->GetYaxis()->GetBinCenter(max_ybin);

  int entry = _houghPhiZHist.at(_houghPhiZHist.size()-1)->GetBinContent(max_xbin,max_ybin);
  _houghPhiZHist.at(_houghPhiZHist.size()-1)->GetXaxis()->SetRange(0,_parHoughNStepTheta);
  if( entry < _parCutHoughPhiZPeak ){
    log("debug") << " => not identified as maximum point in "
		<< _parCutTrackRangeThetaMin << " <theta<" << _parCutTrackRangeThetaMax
		<< " : rho = "    << rho << ", theta = " << theta
		<< ", N_entry = " << entry
		<< std::endl;
    _houghPhiZHist.at(_houghPhiZHist.size()-1)->GetXaxis()->SetRangeUser(0., 180.);
    _houghPhiZHist.at(_houghPhiZHist.size()-1)->GetMaximumBin(max_xbin, max_ybin, max_zbin);
    theta = _houghPhiZHist.at(_houghPhiZHist.size()-1)->GetXaxis()->GetBinLowEdge(max_xbin);
    rho   = _houghPhiZHist.at(_houghPhiZHist.size()-1)->GetYaxis()->GetBinCenter(max_ybin);
  }

  double infinity = 1.0e+10;
  double par1 = ( theta >= 3.0 ? -1.0/TMath::Tan(theta*TMath::Pi()/180.0)*180.0/TMath::Pi() : infinity                        );
  double par0 = ( theta >= 3.0 ?  rho/TMath::Sin(theta*TMath::Pi()/180.0)                   : -rho*infinity*TMath::Pi()/180.0 );

  entry = _houghPhiZHist.at(_houghPhiZHist.size()-1)->GetBinContent(max_xbin,max_ybin);

  log("debug") << "next line-ID is "
	      << _houghPhiZLineNoMax+1 << " : rho(max) = "
	      << rho                      << ", theta = "
	      << theta                    << ", par0 = "
	      << par0                     << ", par1 = "
	      << par1
	      << std::endl;
  
  // break point;  
  if( entry < _parCutHoughPhiZPeak // condition for normal tracks
      || ((theta < _parCutTrackRangeThetaMin || theta > _parCutTrackRangeThetaMax) && entry < _parCutHoughPhiZPeakDeltaRay) // condition for deltaray
      ){
    log("debug") << " => not identified as maximum point : N_entry = " << entry << std::endl;
    return -1;
  }
  //+++++++++++++++++++++++++++++++
  log("debug") << " => identified as maximum point : N_entry = " << entry << std::endl;
  int lineNo = ++_houghPhiZLineNoMax;

  TF1* func_hough_phiz;
  if( TMath::Abs(par0) > 1.0e+5 ) func_hough_phiz = new TF1( Form("func_hough_phiz_%d",lineNo), "[0]+[1]*x", -par0/par1-1.0e-6, -par0/par1+1.0e-6 ); // deltaray
  else                            func_hough_phiz = new TF1( Form("func_hough_phiz_%d",lineNo), "[0]+[1]*x",               0.0,     2*TMath::Pi() ); // normal tracks

  func_hough_phiz->SetLineWidth( 1 );
  func_hough_phiz->SetLineStyle( 1 );
  func_hough_phiz->SetParameter( 0, par0);
  func_hough_phiz->SetParameter( 1, par1);
  _houghPhiZFunc[lineNo] = func_hough_phiz;

  _houghPhiZHist.at(_houghPhiZHist.size()-1)->SetTitle( Form("%s, (#theta,#rho) = (%.1f, %.1f)",_houghPhiZHist.at(_houghPhiZHist.size()-1)->GetTitle(),theta,rho) );

  return 1;
}

int TrkFinding::calcHoughResidualPhiZ(){
  log("debug") << "*****[START : calcHoughResidualPhiZ]*****" << std::endl;
  //if( !_houghPhiZFunc.size() ) return -1; // check if new straight line is found in the latest Hough transformation

  int lineNo = _houghPhiZLineNoMax;
  TH1D* hist_resi_phiz = new TH1D( Form("hist_resi_phiz_%d",lineNo),  Form("Residual_%d(#phi-Z);",lineNo),  100, -20, 20 );
  double par0 = _houghPhiZFunc[lineNo]->GetParameter(0);
  double par1 = _houghPhiZFunc[lineNo]->GetParameter(1);
  log("debug") << "lineNo = " << lineNo << " : par0 = " << par0 << ", par1 = " << par1 << std::endl;

  int cnt_seed = 0;
  double residual;

  for( int ivec=0; ivec<_hitInfoRecoHits.size(); ivec++ ){ // Begin of Loop for Hits
    if( _timeWindowMin > _hitInfoRecoHits.at(ivec)->_time || _timeWindowMax < _hitInfoRecoHits.at(ivec)->_time ) continue; // time-window is "min <= time <= max"
    residual = ( par1 < 1.e+5 ? _hitInfoRecoHits.at(ivec)->_pos.Z() - ( par0+par1*_hitInfoPhi.at(ivec) ) : _hitInfoPhi.at(ivec) - (-par0/par1) );
    
    hist_resi_phiz->Fill( residual );
    if( _hitInfoClusterNo.at(ivec) >= -1 ) continue; // to be checked : -1 or 0
    //if( _hitInfoRecoHits.at(ivec)->_time < _timeWindowMin || _hitInfoRecoHits.at(ivec)->_time > _timeWindowMax ) continue;
    if( (TMath::Abs(residual) < _parCutHoughPhiZSeedResi && par1 < 1.e+5) || (TMath::Abs(residual) < _parCutHoughPhiZSeedResiDeltaRay && par1 >= 1.e+5) ){ // normal-track and deltaray-track
      _hitInfoLineNoPhiZ.at(ivec) = lineNo;
      cnt_seed++;
      log("debug") << Form( "   index = %2d, vaneID = %2d, ",  ivec, _hitInfoVaneID.at(ivec)    )
		   << Form( "(X,Y,Z) = (%7.2f,%7.2f,%7.2f), ", _hitInfoRecoHits.at(ivec)->_pos.X(),_hitInfoRecoHits.at(ivec)->_pos.Y(),_hitInfoRecoHits.at(ivec)->_pos.Z())
		   << Form( "(R,phi) = (%6.2f,%.2f)",          _hitInfoRecoHits.at(ivec)->_pos.Perp(),_hitInfoPhi.at(ivec)            )
		   << std::endl;
    }
  } // End of Loop for Hits
  log("debug") << " => " << cnt_seed << " points are found as seeds" << std::endl;

  _houghPhiZHistResi.push_back( hist_resi_phiz );

  return 1;
}

int TrkFinding::clusteringShort(){
  log("debug") << "*****[START: clusteringShort]*****" << std::endl;
  int lineNo = _houghPhiZLineNoMax;
  if( _houghPhiZFunc[lineNo]->GetParameter(1) < 1.0e+5 ) return clustering3DShort(); //  normal -track
  //if( _houghPhiZFunc[lineNo]->GetParameter(1) < 1.0e+5 ) return clustering3DShortOld(); //  normal -track
  else return -1;
}

///*
int TrkFinding::clustering3DShortOld(){
  log("debug") << "*****[START : Clustering for normal track]*****" << std::endl;
  int lineNo = _houghPhiZLineNoMax;
  bool fl_new_cluster = !true;
  log("debug") << "lineNo = " << lineNo << std::endl;
  
  std::vector<int> cluster_index;
  int index;
  for( unsigned int ivane=0; ivane<_orderVaneID.size(); ivane++ ){
    index = _orderVaneID.at(ivane);
    if( _hitInfoLineNoPhiZ.at(index)!=lineNo ) continue; // select hits associated specific straight line
    if( _hitInfoClusterNo.at(index)       >=0      ) continue;
    if( _hitInfoRecoHits.at(index)->_time < _timeWindowMin || _hitInfoRecoHits.at(index)->_time > _timeWindowMax ) continue; // select hits within specific time window
    _hitInfoClusterNo.at( index ) = -1;
    cluster_index.push_back(index);
    log("debug") << "vane-ID = " << _hitInfoVaneID.at(index) << std::endl;
  }
  log("debug") << "seed cluster : " << cluster_index.size() << std::endl;
  if( cluster_index.size()<3 ){
    log("debug") << " => insufficient seed" << std::endl;
    return -1;
  }
  
  // forward ++++++++++++++++++++++++++;
  double x0;
  double y0;
  double r;
  double dphi;
  double extrap_r;
  double extrap_z;
  std::vector<int> extrap_index;
  int index1, index2, index3;
  int cnt_miss_forward, cnt_miss_backward;
  int cnt_nocrosspoint;
  int pre_target_VaneID;
  int current_VaneID, target_VaneID, result;
  double dev_min;
  int index_min;
  bool match;
  double dev_r, dev_z;
  int cnt_seq;

  for(unsigned int i=0; i<cluster_index.size()-2; i++){
    if( _hitInfoClusterNo.at(cluster_index.at(i))>=0 ) continue;
    for(unsigned int j=i+1; j<cluster_index.size()-1; j++){
      if( _hitInfoClusterNo.at(cluster_index.at(j))>=0 ) continue;
      for(unsigned int k=j+1; k<cluster_index.size(); k++){
	if( _hitInfoClusterNo.at(cluster_index.at(k))>=0 ) continue;
	index1 = cluster_index.at(i);
	index2 = cluster_index.at(j);
      	index3 = cluster_index.at(k);

	if(_hitInfoVaneID.at(index1)==_hitInfoVaneID.at(index2) || _hitInfoVaneID.at(index2)==_hitInfoVaneID.at(index3)) continue;
	if( (_hitInfoVaneID.at(index2)-_hitInfoVaneID.at(index1)+_parGeomNVane)%_parGeomNVane>1 ||
	    (_hitInfoVaneID.at(index3)-_hitInfoVaneID.at(index2)+_parGeomNVane)%_parGeomNVane>1 ) continue;

	extrap_index.clear();
	extrap_index.push_back(index1);
	extrap_index.push_back(index2);
	extrap_index.push_back(index3);

	cnt_miss_forward  = 0;
	cnt_nocrosspoint  = 0;
	pre_target_VaneID = -999;
	log("debug") << "======== start extrapolation ========" << std::endl;

	log("debug") << "******** start extrapolation along forward direction *******" << std::endl;
	while( cnt_miss_forward < _parCutExtrapMiss ){
	  current_VaneID = _hitInfoVaneID.at(index3);
	  target_VaneID  = (pre_target_VaneID < 0 ? _hitInfoVaneID.at(index3)+1 : pre_target_VaneID+1 );
	  log("debug") << "      Extrapolate(forward) from "
		       << "vaneID=" << _hitInfoVaneID.at( index1 ) << "(index=" << index1 << ") & "
		       << "vaneID=" << _hitInfoVaneID.at( index2 ) << "(index=" << index2 << ") & "
		       << "vaneID=" << _hitInfoVaneID.at( index3 ) << "(index=" << index3 << ")" << std::endl;
	  
	  // Find next hit point by extrapolation
	  while(1){
	    if( target_VaneID==_parGeomNVane ) target_VaneID = 0;
	    if( current_VaneID==target_VaneID ){
	      target_VaneID = -999;
	      break;
	    }else{
	      result = extrapolation( index1, index2, index3, target_VaneID, extrap_r, extrap_z, x0, y0, r, dphi );
	      if(result==2){ // out of z
		target_VaneID = -999; // no further extrapolation in this direction
		log("debug") << " => no further extrapolation due to out of z" << std::endl;
		break;
	      }else if(result==1){ // cross point exist
		break;
	      }else{ // no cross point
		cnt_nocrosspoint++;
		target_VaneID++;
		if( cnt_nocrosspoint>=_parCutExtrapNoCross ){
		  target_VaneID = -999;
		  log("debug") << " => no cross point => break" << std::endl;
		  break;
		}
		log("debug") << " => no cross point => continue" << std::endl;
	      }
	    }
	  }
	  if( target_VaneID < 0 ) break;
	  pre_target_VaneID = target_VaneID;
	  log("debug") << "         target VaneID = "
		      << target_VaneID << " : R = "
		      << extrap_r      << ", Z = "
		      << extrap_z      << " : x0 = "
		      << x0            << ", y0 = "
		      << y0            << ", r = "
		      << r             << ", dphi = "
		      << dphi          << std::endl;
	  
	  // compare actual hits with extrapolated point.
	  dev_min   = 10000;
	  index_min = -999;
	  match = false;
	  for( int ivane=0; ivane<_orderVaneID.size(); ivane++ ){
	    if( _hitInfoVaneID.at(_orderVaneID.at(ivane))!=target_VaneID ) continue;
	    if( _hitInfoClusterNo.at(_orderVaneID.at(ivane))>=0          ) continue;
	    for(unsigned int l=0; l<extrap_index.size(); l++){
	      if( extrap_index.at(l)==_orderVaneID.at(ivane) ){
		match = true;
		break;
	      }
	    }
	    if(match) continue;
	    dev_r = TMath::Abs( extrap_r - _hitInfoRecoHits.at(_orderVaneID.at(ivane))->_pos.Perp() );
	    dev_z = TMath::Abs( extrap_z - _hitInfoRecoHits.at(_orderVaneID.at(ivane))->_pos.Z()    );
	    log("debug") << "            index = " << _orderVaneID.at(ivane) << " dev(r) = " << dev_r << ", dev(z) = " << dev_z << ", dev_min = " << sqrt(pow(dev_r,2)+pow(dev_z,2)) << std::endl;
	    if( dev_min > sqrt(pow(dev_r,2)+pow(dev_z,2)) ){
	      dev_min   = sqrt(pow(dev_r,2)+pow(dev_z,2));
	      index_min = _orderVaneID.at(ivane);
	    }
	  }
	  if( dev_min < _parCutExtrapToleranceCoeffDphi*TMath::Abs(dphi) || dev_min < _parCutExtrapTolerance ){
	    extrap_index.push_back(index_min);
	    index1 = index2;
	    index2 = index3;
	    index3 = index_min;
	    cnt_miss_forward = 0;
	    log("debug") << "         => added the hit into the cluster" << std::endl;
	  }else{
	    if( !(extrap_r < _geomRInner || extrap_r > _geomROuter) ) cnt_miss_forward++;
	    log("debug") << "         => can not find hit points by extrapolation : dev_min = " << dev_min << ", cnt_miss_forward = " << cnt_miss_forward << std::endl;
	  }
	}
      
	// backward ++++++++++++++++++++++++++;
	pre_target_VaneID = -999;
	index1 = cluster_index.at(k);
	index2 = cluster_index.at(j);
	index3 = cluster_index.at(i);

	cnt_miss_backward = 0;
	cnt_nocrosspoint = 0;

	log("debug") << "******** start extrapolation along backward direction *******" << std::endl;
 	while( cnt_miss_backward < _parCutExtrapMiss ){
	  current_VaneID = _hitInfoVaneID.at(index3);
	  target_VaneID  = (pre_target_VaneID < 0 ? _hitInfoVaneID.at(index3)-1 : pre_target_VaneID-1 );
	  log("debug") << "      Extrapolate(backward) from "
		       << "vaneID=" << _hitInfoVaneID.at( index1 ) << "(index=" << index1 << ") & "
		       << "vaneID=" << _hitInfoVaneID.at( index2 ) << "(index=" << index2 << ") & "
		       << "vaneID=" << _hitInfoVaneID.at( index3 ) << "(index=" << index3 << ")" << std::endl;
	  
	  // Find next hit point by extrapolation
	  while(1){
	    if( target_VaneID==-1 ) target_VaneID = _parGeomNVane-1;
	    if( current_VaneID==target_VaneID ){
	      target_VaneID = -999;
	      break;
	    }else{
	      result = extrapolation( index1, index2, index3, target_VaneID, extrap_r, extrap_z, x0, y0, r, dphi );
	      if(result==2){ // out of z
		target_VaneID = -999; // no further extrapolation in this direction
		log("debug") << " => no further extrapolation due to out of z" << std::endl;
		break;
	      }else if(result==1){ // cross point exist
		break;
	      }else{ // no cross point
		target_VaneID--;
		cnt_nocrosspoint++;
		if( cnt_nocrosspoint>=_parCutExtrapNoCross ){
		  target_VaneID = -999;
		  log("debug") << " => no cross point => break" << std::endl;
		  break;
		}
		log("debug") << " => no cross point => continue" << std::endl;
	      }
	    }
	  }
	  if( target_VaneID < 0 ) break;
	  pre_target_VaneID = target_VaneID;
	  log("debug") << "         target VaneID = "
		       << target_VaneID << " : R = "
		       << extrap_r      << ", Z = "
		       << extrap_z      << " : x0 = "
		       << x0            << ", y0 = "
		       << y0            << ", r = "
		       << r             << ", dphi = "
		       << dphi          << std::endl;
	  
	  // compare actual hits with extrapolated point.
	  dev_min   = 10000;
	  index_min = -999;
	  match = false;
	  for( int ivane=0; ivane<_orderVaneID.size(); ivane++ ){
	    if( _hitInfoVaneID.at(_orderVaneID.at(ivane))!=target_VaneID ) continue;
	    if( _hitInfoClusterNo.at(_orderVaneID.at(ivane))>=0         ) continue;
	    for(unsigned int l=0; l<extrap_index.size(); l++){
	      if( extrap_index.at(l)==_orderVaneID.at(ivane) ){
		match = true;
		break;
	      }
	    }
	    if( match ) continue;
	    dev_r = TMath::Abs( extrap_r - _hitInfoRecoHits.at(_orderVaneID.at(ivane))->_pos.Perp() );
	    dev_z = TMath::Abs( extrap_z - _hitInfoRecoHits.at(_orderVaneID.at(ivane))->_pos.Z()    );
	    log("debug") << "            index = " << _orderVaneID.at(ivane) << " dev(r) = " << dev_r << ", dev(z) = " << dev_z << ", dev_min = " << sqrt(pow(dev_r,2)+pow(dev_z,2)) << std::endl;
	    if( dev_min > sqrt(pow(dev_r,2)+pow(dev_z,2)) ){
	      dev_min   = sqrt(pow(dev_r,2)+pow(dev_z,2));
	      index_min = _orderVaneID.at(ivane);
	    }
	  }
	  if( dev_min < _parCutExtrapToleranceCoeffDphi*TMath::Abs(dphi) || dev_min < _parCutExtrapTolerance ){

	    extrap_index.insert(extrap_index.begin(), index_min);
	    cnt_miss_backward = 0;
	    index1 = index2;
	    index2 = index3;
	    index3 = index_min;
	    log("debug") << "         => added the hit into the cluster" << std::endl;
	  }else{
	    if( !(extrap_r < _geomRInner || extrap_r > _geomROuter) ) cnt_miss_backward++;
	    log("debug") << "         => can not find hit points by extrapolation : dev_min = " << dev_min << ", cnt_miss_backward = " << cnt_miss_backward << std::endl;
	  }
	}

	if( extrap_index.size()>4 ){
	  _clusterNoMax++;
	  fl_new_cluster = true;
	  log("debug") << " ===== clusterNo = " << _clusterNoMax << " =====" << std::endl;
	  for( unsigned int l=0; l<extrap_index.size(); l++ ){
	    log("debug") << "index = "       << extrap_index.at(l)
			<< ", vaneID = "    << _hitInfoVaneID.at(extrap_index.at(l))
			<< ", clusterNo = " << _clusterNoMax
			<< std::endl;
	    _hitInfoClusterNo.at( extrap_index.at(l) ) = _clusterNoMax;
	  }
	  // assign sequence number in the cluster
	  cnt_seq = 0;
	  std::vector<Int_t> tmp_order_cluster;
	  log("debug") << "       ==> clusterNo = " << _clusterNoMax << ", lineNo = " << lineNo << " with " << extrap_index.size() << " hits" << std::endl;
	  if( TMath::Abs(_hitInfoRecoHits.at(_orderVaneID.at(extrap_index.at(0)))->_pos.Z()) < TMath::Abs(_hitInfoRecoHits.at(_orderVaneID.at(extrap_index.at(extrap_index.size()-1)))->_pos.Z()) ){
	    for( int ihit=0; ihit<extrap_index.size(); ihit++ ){
	      _hitInfoSequenceNo.at( extrap_index.at(ihit) ) = cnt_seq;
	      tmp_order_cluster.push_back(extrap_index.at(ihit));
	      cnt_seq++;
	    }
	  }else{
	    for( int ihit=extrap_index.size()-1; ihit>=0; ihit-- ){
	      _hitInfoSequenceNo.at( extrap_index.at(ihit) ) = cnt_seq;
	      tmp_order_cluster.push_back(extrap_index.at(ihit));
	      cnt_seq++;
	    }
	  }
	  _mapClusteringResult      [_clusterNoMax] = tmp_order_cluster;
	  _mapActiveClusterNo       [_clusterNoMax] = _timeWindowMax;
	  _mapActiveClusterNoLineNo[_clusterNoMax] = lineNo;
	  _setActiveLineNo.insert(lineNo);
	  break;
	}else{
	  log("debug") << " => not identified as cluster due to small number of hits" << std::endl;
	}
      } // k
      if( _hitInfoClusterNo.at(cluster_index.at(i))>=0 ) break;
    } // j
  } // i

  if( !fl_new_cluster ){
    log("debug") << "       ==> New cluster is not found from lineNo = " << lineNo << std::endl;
  }

  log("debug") << "   *****[FINISH : Clustering]*****" << std::endl;

  return 1;
}
//*/

int TrkFinding::clustering3DShort(){
  log("debug") << "*****[START : Clustering for normal track]*****" << std::endl;
  int lineNo = _houghPhiZLineNoMax;
  bool fl_new_cluster = !true;
  log("debug") << "lineNo = " << lineNo << std::endl;
  
  std::vector<int> cluster_index; // seed hits (hit close to hough-fit line)
  //int index;
  //for( unsigned int ivane=0; ivane<_orderVaneID.size(); ivane++ ){
  //  index = _orderVaneID.at(ivane);
  for(int index=0; index<getNhits(); index++){
    if( _hitInfoLineNoPhiZ.at(index)!=lineNo ) continue; // select hits associated specific straight line
    if( _hitInfoClusterNo.at(index) >=0      ) continue;
    if( _timeWindowMin > _hitInfoRecoHits.at(index)->_time || _timeWindowMax < _hitInfoRecoHits.at(index)->_time ) continue; // time-window is "min <= time <= max"
    _hitInfoClusterNo.at( index ) = -1;
    cluster_index.push_back(index);
    log("debug") << "vane-ID = " << _hitInfoVaneID.at(index) << std::endl;
  }
  log("debug") << "seed cluster : " << cluster_index.size() << std::endl;
  if( cluster_index.size()<3 ){
    log("debug") << " => insufficient seed" << std::endl;
    return -1;
  }
  
  std::vector<int> extrap_index;
  int index1, index2, index3;
  int cnt_seq;

  for(unsigned int i=0; i<cluster_index.size()-2; i++){
    index1 = cluster_index.at(i);
    if( _hitInfoClusterNo.at(index1)>=0 ) continue;
    for(unsigned int j=i+1; j<cluster_index.size()-1; j++){
      index2 = cluster_index.at(j);
      if( _hitInfoClusterNo.at(index2)>=0 ) continue;
      if(_hitInfoVaneID.at(index1)==_hitInfoVaneID.at(index2) ) continue;
      if( (_hitInfoVaneID.at(index2)-_hitInfoVaneID.at(index1)+_parGeomNVane)%_parGeomNVane>1 ) continue;
      for(unsigned int k=j+1; k<cluster_index.size(); k++){
	//index1 = cluster_index.at(i);
	//index2 = cluster_index.at(j);
      	index3 = cluster_index.at(k);
	if( _hitInfoClusterNo.at(index3)>=0 ) continue;

	if( _hitInfoVaneID.at(index2)==_hitInfoVaneID.at(index3)) continue;
	if( (_hitInfoVaneID.at(index3)-_hitInfoVaneID.at(index2)+_parGeomNVane)%_parGeomNVane>1 ) continue;

	//if(_hitInfoVaneID.at(index1)==_hitInfoVaneID.at(index2) || _hitInfoVaneID.at(index2)==_hitInfoVaneID.at(index3)) continue;
	//if( (_hitInfoVaneID.at(index2)-_hitInfoVaneID.at(index1)+_parGeomNVane)%_parGeomNVane>1 ||
	//    (_hitInfoVaneID.at(index3)-_hitInfoVaneID.at(index2)+_parGeomNVane)%_parGeomNVane>1 ) continue;

	extrap_index.clear();
	extrap_index.push_back(index1);
	extrap_index.push_back(index2);
	extrap_index.push_back(index3);

	expandClustering( extrap_index, true  ); // forward
	expandClustering( extrap_index, false ); // backward

	if( extrap_index.size()>4 ){
	  _clusterNoMax++;
	  fl_new_cluster = true;
	  log("debug") << " ===== clusterNo = " << _clusterNoMax << " =====" << std::endl;
	  for( unsigned int l=0; l<extrap_index.size(); l++ ){
	    log("debug") << "index = "       << extrap_index.at(l)
			 << ", vaneID = "    << _hitInfoVaneID.at(extrap_index.at(l))
			 << ", clusterNo = " << _clusterNoMax
			 << std::endl;
	    _hitInfoClusterNo.at( extrap_index.at(l) ) = _clusterNoMax;
	  }

	  // assign sequence number in the cluster
	  cnt_seq = 0;
	  std::vector<Int_t> tmp_order_cluster;
	  log("debug") << "       ==> clusterNo = " << _clusterNoMax << ", lineNo = " << lineNo << " with " << extrap_index.size() << " hits" << std::endl;
	  //if( TMath::Abs(_hitInfoRecoHits.at(_orderVaneID.at(extrap_index.at(0)))->_pos.Z()) < TMath::Abs(_hitInfoRecoHits.at(_orderVaneID.at(extrap_index.at(extrap_index.size()-1)))->_pos.Z()) ){ // judge which is head/tail of track based on Z position
	  // judge head/tail of track assuming positive charge
	  bool poscharge = true;
	  double x0,y0,r;
	  /*
	  if( circleBy3Point( _hitInfoRecoHits.at(_orderVaneID.at(extrap_index.at(0)))->_pos.X(),
			      _hitInfoRecoHits.at(_orderVaneID.at(extrap_index.at(0)))->_pos.Y(),
			      _hitInfoRecoHits.at(_orderVaneID.at(extrap_index.at(1)))->_pos.X(),
			      _hitInfoRecoHits.at(_orderVaneID.at(extrap_index.at(1)))->_pos.Y(),
			      _hitInfoRecoHits.at(_orderVaneID.at(extrap_index.at(2)))->_pos.X(),
			      _hitInfoRecoHits.at(_orderVaneID.at(extrap_index.at(2)))->_pos.Y(),
			      x0,y0,r)>0 ){
	    double phi = atan2(_hitInfoRecoHits.at(_orderVaneID.at(extrap_index.at(0)))->_pos.Y()-y0,
				_hitInfoRecoHits.at(_orderVaneID.at(extrap_index.at(0)))->_pos.X()-x0);
	    double vx = sin(phi);
	    double vy = -cos(phi);
	    if( ((_hitInfoRecoHits.at(_orderVaneID.at(extrap_index.at(0)))->_pos.X()-x0)*vy
		 -(_hitInfoRecoHits.at(_orderVaneID.at(extrap_index.at(0)))->_pos.Y()-y0)*vx) > 0 ){
	      poscharge = false;
	    }
	  */
	  if( circleBy3Point( _hitInfoRecoHits.at(extrap_index.at(0))->_pos.X(),
			      _hitInfoRecoHits.at(extrap_index.at(0))->_pos.Y(),
			      _hitInfoRecoHits.at(extrap_index.at(1))->_pos.X(),
			      _hitInfoRecoHits.at(extrap_index.at(1))->_pos.Y(),
			      _hitInfoRecoHits.at(extrap_index.at(2))->_pos.X(),
			      _hitInfoRecoHits.at(extrap_index.at(2))->_pos.Y(),
			      x0,y0,r)>0 ){
	    double phi = atan2(_hitInfoRecoHits.at(extrap_index.at(0))->_pos.Y()-y0,
				_hitInfoRecoHits.at(extrap_index.at(0))->_pos.X()-x0);
	    double vx = sin(phi);
	    double vy = -cos(phi);
	    if( ((_hitInfoRecoHits.at(extrap_index.at(0))->_pos.X()-x0)*vy
		 -(_hitInfoRecoHits.at(extrap_index.at(0))->_pos.Y()-y0)*vx) > 0 ){
	      poscharge = false;
	    }
	  }
	  if( poscharge ){
	    for( int ihit=0; ihit<extrap_index.size(); ihit++ ){
	      _hitInfoSequenceNo.at( extrap_index.at(ihit) ) = cnt_seq;
	      tmp_order_cluster.push_back(extrap_index.at(ihit));
	      cnt_seq++;
	    }
	  }else{
	    for( int ihit=extrap_index.size()-1; ihit>=0; ihit-- ){
	      _hitInfoSequenceNo.at( extrap_index.at(ihit) ) = cnt_seq;
	      tmp_order_cluster.push_back(extrap_index.at(ihit));
	      cnt_seq++;
	    }
	  }
	  _mapClusteringResult      [_clusterNoMax] = tmp_order_cluster;
	  _mapActiveClusterNo       [_clusterNoMax] = _timeWindowMax;
	  _mapActiveClusterNoLineNo[_clusterNoMax] = lineNo;
	  _setActiveLineNo.insert(lineNo);
	  break;
	}else{
	  log("debug") << " => not identified as cluster due to small number of hits" << std::endl;
	}
      } // k
      if( _hitInfoClusterNo.at(cluster_index.at(i))>=0 ) break;
    } // j
  } // i

  if( !fl_new_cluster ) log("debug") << "       ==> New cluster is not found from lineNo = " << lineNo << std::endl;
  log("debug") << "   *****[FINISH : Clustering]*****" << std::endl;

  return 1;
}


void TrkFinding::judgeFinding( std::vector<std::pair<int,int> > *result ){
  log("debug") << "*****[START: judgeFinding]*****" << std::endl;
  log("debug") << "clusterNo_max = " << _clusterNoMax << std::endl;
  if(_clusterNoMax<0) return;

  std::vector<std::vector<int> > clusterNo_index;
  for( unsigned int i=0; i<=_clusterNoMax; i++ ){
    std::vector<int> v;
    clusterNo_index.push_back(v);
  }
  
  // fill index by T order
  int index;
  int clusterNo;
  for( unsigned int i=0; i<_orderRecTime.size(); i++ ){ // better to use _orderMCTruthTime ?
    index     = _orderRecTime.at(i);
    clusterNo = _hitInfoClusterNo.at(index);
    if( _hitInfoClusterNo.at(index) >= 0 ) clusterNo_index.at( clusterNo ).push_back( index );
  }

  int  preVaneID  = -999;
  int  preMuonID  = -999; 
  bool increasing = kFALSE;
  bool decreasing = kFALSE;
  int  count      = 0;
  int  countMax   = 0;
  int  MuonIDMax  = 0;

  // search for continous hits
  for( unsigned int i=0; i<clusterNo_index.size(); i++ ){ // BEGIN LOOP for clusterNo
    log("debug") << "clusterNo = " << i << "********************************" << std::endl;
    preVaneID = -999;
    preMuonID = -999;
    increasing = kFALSE;
    decreasing = kFALSE;
    count      = 0;
    countMax   = 0;
    MuonIDMax  = -1;
    for( unsigned int j=0; j<clusterNo_index.at(i).size(); j++ ){ // BEGIN LOOP for sequenceNo
      index = clusterNo_index.at(i).at(j);
      log("debug") << "      index = " << index << ", vane-ID = " << _hitInfoVaneID.at(index) << " : ";
      if( _hitInfoVaneID.at(index)==preVaneID ){
	log("debug") << " => do nothing" << std::endl;
	; // do nothing
      }else if( increasing &&
	       _hitInfoVaneID.at(index)==((preVaneID+1)%_parGeomNVane) &&
	       _hitInfoMuonID.at(index)==preMuonID ){
	log("debug") << " => increasing vane-ID" << std::endl;
	preVaneID = _hitInfoVaneID.at(index);
	count++;
	decreasing = kFALSE;
      }else if( decreasing && 
		_hitInfoVaneID.at(index)==((preVaneID+_parGeomNVane-1)%_parGeomNVane) &&
		_hitInfoMuonID.at(index)==preMuonID ){
	log("debug") << " => decreasing vane-ID" << std::endl;
	preVaneID = _hitInfoVaneID.at(index);
	count++;
	increasing = kFALSE;
      }else{
	log("debug") << " => missing hit or different muon origin hit" << std::endl;
	preVaneID = _hitInfoVaneID.at(index);
	preMuonID = _hitInfoMuonID.at(index);
	count = 1;
	increasing = kTRUE;
	decreasing = kTRUE;
      }

      if( count>countMax && preMuonID >=0 ){
	countMax  = count;
	MuonIDMax = preMuonID;
      }
    } // END LOOP for sequenceNo
    
    if(countMax>=_parThresholdSuccess){
      result->push_back(std::make_pair((int)i, MuonIDMax));
      log("debug") << "clusterNo:" << i << " finding is succeeded for ID:" << MuonIDMax << std::endl; 
    }
  } // END LOOP for clusterNo

}


void TrkFinding::evaluateTrkFinding(){
  log("debug") << "*****[START: evaluateTrkFinding]*****" << std::endl;

  _performance_nhits.clear();
  _performance_ngoodhits.clear();
  _performance_ncontinuous_hits.clear();
  _performance_muonid.clear();
  _performance_clusterno.clear();
  _performance_existmuonid.clear();
  _performance_existenergy.clear();
  _performance_existpt.clear();
  _performance_existpz.clear();

  std::map<int, TGraph*> map_g_mc_hit_xy;    // key is muon-ID
  std::map<int, TGraph*> map_g_rec_hit_xy;   // key is clusterNo
  std::map<int, TGraph*> map_g_mc_hit_phiz;  // key is muon-ID
  std::map<int, TGraph*> map_g_rec_hit_phiz; // key is clusterNo

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // check MC-truth information
  std::map<int,const MCParticle*>::iterator it_muon     = _mapMCMuon.begin();
  std::map<int,const MCParticle*>::iterator it_positron = _mapMCPositron.begin();
  log("info") << _mapMCMuon.size() << " muons exist in this event" << std::endl;

  while(it_muon != _mapMCMuon.end() ||
	it_positron != _mapMCPositron.end()
	){
    const MCParticle* mc_muon     = it_muon->second;
    const MCParticle* mc_positron = it_positron->second;
    int    muon_id                = it_muon->first;
    double energy                 = mc_positron->_p.E();
    log("debug") << "    Muon-ID = " << muon_id
		<< " : Energy = "   << energy << " MeV"
      		<< ", Pt = "        << mc_positron->_p.Perp() << " MeV"
      		<< ", Pz = "        << mc_positron->_p.Z()    << " MeV"
      		<< std::endl;

    _performance_existmuonid.push_back(muon_id);
    _performance_existenergy.push_back(energy);
    _performance_existpt.push_back(mc_positron->_p.Perp());
    _performance_existpz.push_back(mc_positron->_p.Z());

    map_g_mc_hit_xy   [muon_id] = makeGraphMuonIDXY  (muon_id, ( _hitInfoRecoHits.size() ? _hitInfoRecoHits.at(_orderRecTime.at(0))->_time : -999 ), ( _hitInfoRecoHits.size() ? _hitInfoRecoHits.at(_orderRecTime.at(_hitInfoRecoHits.size()-1))->_time : -999 ) );
    map_g_mc_hit_phiz [muon_id] = makeGraphMuonIDPhiZ(muon_id, ( _hitInfoRecoHits.size() ? _hitInfoRecoHits.at(_orderRecTime.at(0))->_time : -999 ), ( _hitInfoRecoHits.size() ? _hitInfoRecoHits.at(_orderRecTime.at(_hitInfoRecoHits.size()-1))->_time : -999 ) );

    ++it_muon;
    ++it_positron;
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
  // check reconstructed cluster
  log("debug") << "clusterNo_max = " << _clusterNoMax << std::endl;
  if(_clusterNoMax<0) return;
  std::vector<int> candidatesPerMuon(_mapMCMuon.size(),0);

  std::map<int,std::vector<int> >::iterator it_cluster = _mapClusteringResult.begin();
  std::multiset<int> cnt_muonid;    // muonid
  std::set<int>      index_rec;     // index
  std::vector<int>   index_mctruth; // index

  while( it_cluster != _mapClusteringResult.end() ){ // BEGIN LOOP for clusterNo
    int clusterNo = it_cluster->first;
    map_g_rec_hit_xy  [clusterNo] = makeGraphClusterNoXY  ( clusterNo, ( _hitInfoRecoHits.size() ? _hitInfoRecoHits.at(_orderRecTime.at(0))->_time : -999 ), ( _hitInfoRecoHits.size() ? _hitInfoRecoHits.at(_orderRecTime.at(_hitInfoRecoHits.size()-1))->_time : -999 ) );
    map_g_rec_hit_phiz[clusterNo] = makeGraphClusterNoPhiZ( clusterNo, ( _hitInfoRecoHits.size() ? _hitInfoRecoHits.at(_orderRecTime.at(0))->_time : -999 ), ( _hitInfoRecoHits.size() ? _hitInfoRecoHits.at(_orderRecTime.at(_hitInfoRecoHits.size()-1))->_time : -999 ) );

    std::vector<int> sequence = it_cluster->second;
    cnt_muonid.clear();
    index_rec.clear();
    for( int iseq=0; iseq<sequence.size(); iseq++ ){ // BEGIN LOOP for index
      int index = sequence.at(iseq);
      log("debug") << "      index = " << setw(4) << right << index
		   << ", rec-time = "  << setw(6) << right << Form("%.2f",_hitInfoRecoHits.at(index)->_time)
		   << ", mc-time = "   << setw(8) << right << Form("%.3f",getMCTruthTimefromRecoHit(_hitInfoRecoHits.at(index)))
		   << ", x = "         << setw(7) << right << Form("%.2f",_hitInfoRecoHits.at(index)->_pos.X())
		   << ", y = "         << setw(7) << right << Form("%.2f",_hitInfoRecoHits.at(index)->_pos.Y())
		   << ", z = "         << setw(7) << right << Form("%.2f",_hitInfoRecoHits.at(index)->_pos.Z())
		   << ", R = "         << setw(7) << right << Form("%.2f",_hitInfoRecoHits.at(index)->_pos.Perp())
		   << ", Phi = "       << setw(5) << right << Form("%.2f",_hitInfoPhi.at(index))
		   << ", vane-ID = "   << setw(2) << right << _hitInfoVaneID.at(index)
		   << ", muon-ID = "   << setw(3) << right << _hitInfoMuonID.at(index)
		   << ", trueType = "      << _hitInfoRecoHits.at(index)->_trueType
		   << std::endl;
      cnt_muonid.insert(_hitInfoMuonID.at(index));
      index_rec.insert(index);
    } // END LOOP for index

    // check which muon-id is dominant.
    int muonid  = getDominantMuonID(cnt_muonid);
    int cnt_max = cnt_muonid.count(muonid);

    double purity = (double)cnt_max/cnt_muonid.size();
    _performance_nhits.push_back(cnt_muonid.size());
    _performance_ngoodhits.push_back(cnt_max);
    _performance_muonid.push_back(muonid);
    _performance_clusterno.push_back(clusterNo);
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // check mc-truth information
    index_mctruth.clear();
    if( muonid >= 0 ){
      log("debug") << "+++++++++++++++++++ MC-truth information for muon-ID = " << muonid
		   << ", X = " << _mapMCPositron.at(muonid)->_prodVertex.X()
		   << ", Y = " << _mapMCPositron.at(muonid)->_prodVertex.Y()
		   << ", Z = " << _mapMCPositron.at(muonid)->_prodVertex.Z()
		   << std::endl;
    }
    int max_continuous_hits = 0;
    int cnt_continuous_hits = 0;
    for( int ivec=0; ivec<_orderMCTruthTime.size(); ivec++ ){
      int index = _orderMCTruthTime.at(ivec);
      if( _hitInfoRecoHits.at(index)->_trueType ) continue; // remove ghost or noise hits
      if( _hitInfoMuonID.at(index)!=muonid  ) continue;
      index_mctruth.push_back(index);
      log("debug") << "   index = "    << setw(4) << right << index
		   << ", vane-ID = "   << setw(2) << right << _hitInfoVaneID.at(index)
		   << ", rec-time = "  << setw(6) << right << Form("%.2f",_hitInfoRecoHits.at(index)->_time)
		   << ", mc-time = "   << setw(8) << right << Form("%.3f",getMCTruthTimefromRecoHit(_hitInfoRecoHits.at(index)))
		   << ", x = "         << setw(7) << right << Form("%.2f",_hitInfoRecoHits.at(index)->_pos.X())
		   << ", y = "         << setw(7) << right << Form("%.2f",_hitInfoRecoHits.at(index)->_pos.Y())
		   << ", z = "         << setw(7) << right << Form("%.2f",_hitInfoRecoHits.at(index)->_pos.Z())
		   << ", R = "         << setw(7) << right << Form("%.2f",_hitInfoRecoHits.at(index)->_pos.Perp())
		   << ", Phi = "       << setw(5) << right << Form("%.2f",_hitInfoPhi.at(index))
		   << " : isdetected = " << index_rec.count(index)
		   << std::endl;
      if( index_rec.count(index) ){
	cnt_continuous_hits++;
      }else{
	if( cnt_continuous_hits > max_continuous_hits ) max_continuous_hits = cnt_continuous_hits;
	cnt_continuous_hits = 0;
      }
      //checkRecoHits(_hitInfoRecoHits.at(index));
    }

    log("info") << "clusterNo = "                  << setw(3) << right << clusterNo
		<< " with "                        << setw(3) << right << it_cluster->second.size() << " hits"
		<< " seems coming from muon-ID = " << setw(3) << right << muonid
		<< ", purity = "                   << setw(8) << Form("%.2f",purity)
		<< ", " << setw(2) << right << max_continuous_hits << " continuous hits are included." << std::endl;
    _performance_ncontinuous_hits.push_back(max_continuous_hits);
    log("info") << "dominant muonid = " << muonid << endl;

    ++it_cluster;
  } // END LOOP for clusterNo

  if( _indepStudy ) _tree->Fill();
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // check one by one on event display
  log("debug") << "Reconstructed Clusters are checked one by one on event display" <<  std::endl;
  it_muon     = _mapMCMuon.begin();
  it_positron = _mapMCPositron.begin();
  while(it_muon != _mapMCMuon.end() ||
	it_positron != _mapMCPositron.end()
	){
    int muonid = it_muon->first;

    _canvas1->cd(3);
    if( muonid>=0 ){
      gPad->DrawFrame(-350,-350,350,350, Form("MuonID#%d, E = %.1f MeV;X [m];Y [mm]",muonid, _mapMCPositron.at(muonid)->_p.E()) );
      for( int ivane=0; ivane<_parGeomNVane; ivane++ ) _detectorVane[ivane]->Draw("same");
      _detectorMuonOrbit->Draw("Lsame");
      _detectorPole ->Draw("Lsame");
      
      _mapDecayPointXY    [muonid]->Draw("Psame");  
      _mapDecayDirectionXY[muonid]->Draw();
      _mapIdealTrajectory [muonid]->Draw("Psame");  
      map_g_mc_hit_xy[muonid]->SetMarkerStyle(20);
      map_g_mc_hit_xy[muonid]->SetMarkerSize(0.5);
      if( map_g_mc_hit_xy[muonid]->GetN() ) map_g_mc_hit_xy [muonid]->Draw("Psame");
    
      _canvas1->cd(6)->Clear();
      _canvas1->cd(6);
      gPad->DrawFrame(0.0,-250,2.0*TMath::Pi(),250, ";#phi [rad];Z [mm]");
      _mapDecayPointPhiZ    [muonid]->Draw("Psame");  
      _mapDecayDirectionPhiZ[muonid]->Draw();
      map_g_mc_hit_phiz[muonid]->SetMarkerStyle(20);
      map_g_mc_hit_phiz[muonid]->SetMarkerSize(0.5);
      if( map_g_mc_hit_phiz [muonid]->GetN() ) map_g_mc_hit_phiz [muonid]->Draw("Psame");
    }

    it_cluster = _mapClusteringResult.begin();
    int cnt_cluster = 0;
    int cnt_col     = 1;
    while( it_cluster != _mapClusteringResult.end() ){ // BEGIN LOOP for clusterNo
      int clusterNo        = it_cluster->first;
      int estimated_muonid = _performance_muonid.at(cnt_cluster);
      if( muonid == estimated_muonid ){
	_canvas1->cd(3);
	map_g_rec_hit_xy[clusterNo]->SetMarkerStyle(24);
	map_g_rec_hit_xy[clusterNo]->SetMarkerSize(0.7);
	map_g_rec_hit_xy[clusterNo]->SetMarkerColor(1+cnt_col);
	if( map_g_rec_hit_xy[clusterNo]->GetN() ) map_g_rec_hit_xy[clusterNo]->Draw("Psame");
	_canvas1->cd(6);
	map_g_rec_hit_phiz[clusterNo]->SetMarkerStyle(24);
	map_g_rec_hit_phiz[clusterNo]->SetMarkerSize(0.7);
	map_g_rec_hit_phiz[clusterNo]->SetMarkerColor(1+cnt_col);
	if( map_g_rec_hit_phiz[clusterNo]->GetN() ) map_g_rec_hit_phiz[clusterNo]->Draw("Psame");
	cnt_col++;
      }
      cnt_cluster++;
      ++it_cluster;
    }
    
    ++it_muon;
    ++it_positron;
    if( _parDrawLevel ) _canvas1->WaitPrimitive();
  }

  // clear
  std::map<int,TGraph*>::iterator it_g_mc_hit_xy = map_g_mc_hit_xy.begin();
  while( it_g_mc_hit_xy != map_g_mc_hit_xy.end() ){
    delete it_g_mc_hit_xy->second;
    map_g_mc_hit_xy.erase(it_g_mc_hit_xy++);
  }
  std::map<int,TGraph*>::iterator it_g_rec_hit_xy = map_g_rec_hit_xy.begin();
  while( it_g_rec_hit_xy != map_g_rec_hit_xy.end() ){
    delete it_g_rec_hit_xy->second;
    map_g_rec_hit_xy.erase(it_g_rec_hit_xy++);
  }
  std::map<int,TGraph*>::iterator it_g_mc_hit_phiz = map_g_mc_hit_phiz.begin();
  while( it_g_mc_hit_phiz != map_g_mc_hit_phiz.end() ){
    delete it_g_mc_hit_phiz->second;
    map_g_mc_hit_phiz.erase(it_g_mc_hit_phiz++);
  }
  std::map<int,TGraph*>::iterator it_g_rec_hit_phiz = map_g_rec_hit_phiz.begin();
  while( it_g_rec_hit_phiz != map_g_rec_hit_phiz.end() ){
    delete it_g_rec_hit_phiz->second;
    map_g_rec_hit_phiz.erase(it_g_rec_hit_phiz++);
  }

  return;
}

void TrkFinding::writeRecoTrks(){
  log("debug") << "*****[START: writeRecoTrks]*****" << std::endl;
  if(_clusterNoMax<0) return;
  const vector<const MCParticle*> &mc = TreeManager::instance()->getBranchVec<MCParticle>(_parMCParticleName.c_str());

  std::map<int,std::vector<int> >::iterator it_cluster = _mapClusteringResult.begin();
  std::multiset<int> cnt_muonid;

  while( it_cluster != _mapClusteringResult.end() ){ // BEGIN LOOP for cluster
    Track* track  = new Track();
    track->_charge = 1; // assumed to be positron
    int clusterNo = it_cluster->first;

    std::vector<int> sequence = it_cluster->second;
    cnt_muonid.clear();
    for( int iseq=0; iseq<sequence.size(); iseq++ ){
      track->_recoHits.push_back(_hitInfoRecoHits.at(sequence.at(iseq)));
      cnt_muonid.insert(_hitInfoMuonID.at(sequence.at(iseq)));
      if( iseq==0 ){
	track->_pos = _hitInfoRecoHits.at(sequence.at(iseq))->_pos;
	track->_time  = _hitInfoRecoHits.at(sequence.at(iseq))->_time;
      }
    }
    
    // calculate initial momentum
    if( track->_recoHits.size()>=3 ){
      double x0,y0,r;
      if( circleBy3Point( track->_recoHits.at(0)->_pos.X(),
			  track->_recoHits.at(0)->_pos.Y(),
			  track->_recoHits.at(1)->_pos.X(),
			  track->_recoHits.at(1)->_pos.Y(),
			  track->_recoHits.at(2)->_pos.X(),
			  track->_recoHits.at(2)->_pos.Y(),
			  x0,y0,r ) ){
	double pt = 3.*r*0.299792458; // MeV/c
	double phi0 = atan2(track->_recoHits.at(0)->_pos.Y()-y0, track->_recoHits.at(0)->_pos.X()-x0);
	double px = pt*sin(phi0);
	double py = -pt*cos(phi0);
	double phi1 = atan2(track->_recoHits.at(1)->_pos.Y()-y0, track->_recoHits.at(1)->_pos.X()-x0);
	double pz = pt*(track->_recoHits.at(0)->_pos.Z()-track->_recoHits.at(1)->_pos.Z())/(r*(phi0-phi1));
	track->_p = TVector3(px,py,pz);
	log("debug") << "px = " << px << " py = " << py << " pz = " << pz << endl;
      }
    }

    // check which muon-id is dominant.
    int muonid  = getDominantMuonID(cnt_muonid);
    int cnt_max = cnt_muonid.count(muonid);

    if( muonid>=0 ){
      //track->_prodVertex = _mapMCPositron.at(muonid)->_prodVertex; // MC-truth information is input tentatively.
      //track->_p          = _mapMCPositron.at(muonid)->_p.Vect();   // MC-truth information is input tentatively.
      //track->_decayTime  = _mapMCPositron.at(muonid)->_time;       // MC-truth information is input tentatively.
      //track->_charge     = 1; // hard code
      track->_mcp        = const_cast<MCParticle *>(_mapMCPositron.at(muonid));
    }
    _tracks.push_back(track);
    ++it_cluster;
  } // END LOOP for cluster

  return;
}

void TrkFinding::drawEvtDisplay(int nEvent){
  if( !_parDrawLevel ) return;
  log("debug") << "*****[START: drawEvtDisplay]*****" << std::endl;
  makeGraphXY  ( _graphHitPointXYAll,   _graphHitPointXYSig,   _graphHitPointXYBkg   );
  makeGraphPhiZ( _graphHitPointPhiZAll, _graphHitPointPhiZSig, _graphHitPointPhiZBkg );

  _canvas1->cd(1)->Clear();
  _canvas1->cd(1);
  gPad->DrawFrame(-350,-350,350,350, Form("EvtNo.%d : Time Window=[%.1f-%.1f];X [m];Y [mm]",nEvent,_timeWindowMin, _timeWindowMax) ); // error message from here.
  // Error in <TGraphPainter::PaintGraph>: illegal number of points (0) ???
  for( int ivane=0; ivane<_parGeomNVane; ivane++ ) _detectorVane[ivane]->Draw("same");
  _detectorMuonOrbit->Draw("Lsame");
  _detectorPole ->Draw("Lsame");

  std::set<int>::iterator it_active_muonid = _setActiveMuonID.begin();
  while( it_active_muonid != _setActiveMuonID.end() ){
    int muonid = *it_active_muonid;
    log("info") << "muonid = " << muonid << endl;
    if( muonid>=0 ){
      _mapDecayPointXY    [muonid]->Draw("Psame");  
      _mapDecayDirectionXY[muonid]->Draw();
      _mapIdealTrajectory [muonid]->Draw("Psame");  
    }
    ++ it_active_muonid;
  }

  _graphHitPointXYBkg->Draw("Psame");
  _graphHitPointXYSig->Draw("Psame");

  _canvas1->cd(2)->Clear();
  _canvas1->cd(2);
  gPad->DrawFrame(0.0,-250,2.0*TMath::Pi(),250, ";#phi [rad];Z [mm]");

  it_active_muonid = _setActiveMuonID.begin();
  while( it_active_muonid != _setActiveMuonID.end() ){
    int muonid = *it_active_muonid;
    if( muonid>= 0 ){
      _mapDecayPointPhiZ    [muonid]->Draw("Psame");  
      _mapDecayDirectionPhiZ[muonid]->Draw();
    }
    ++ it_active_muonid;
  }

  _graphHitPointPhiZBkg->Draw("Psame");
  _graphHitPointPhiZSig->Draw("Psame");

  _canvas1->cd(4)->Clear();
  _canvas1->cd(4);
  gPad->DrawFrame(-350,-350,350,350, ";X [m];Y [mm]" ); // error message comes from here.
  for( int ivane=0; ivane<_parGeomNVane; ivane++ ) _detectorVane[ivane]->Draw("same");
  _detectorMuonOrbit->Draw("Lsame");
  _detectorPole     ->Draw("Lsame");

  it_active_muonid = _setActiveMuonID.begin();
  while( it_active_muonid != _setActiveMuonID.end() ){
    int muonid = *it_active_muonid;
    if( muonid>=0 ){
      _mapDecayPointXY    [muonid]->Draw("Psame");  
      _mapDecayDirectionXY[muonid]->Draw();
      _mapIdealTrajectory [muonid]->Draw("Psame");  
    }
    ++ it_active_muonid;
  }

  _graphHitPointXYAll->Draw("Psame");

  std::map<int,TGraph*>::iterator it_mapClusteredHitXY = _mapClusteredHitXY.begin();
  while( it_mapClusteredHitXY != _mapClusteredHitXY.end() ){
    it_mapClusteredHitXY->second->Draw("Psame");
    ++it_mapClusteredHitXY;
  }

  _canvas1->cd(5)->Clear();
  _canvas1->cd(5);
  gPad->DrawFrame(0.0,-250,2.0*TMath::Pi(),250, ";#phi [rad];Z [mm]");

  it_active_muonid = _setActiveMuonID.begin();
  while( it_active_muonid != _setActiveMuonID.end() ){
    int muonid = *it_active_muonid;
    if( muonid>=0 ){
      _mapDecayPointPhiZ    [muonid]->Draw("Psame");  
      _mapDecayDirectionPhiZ[muonid]->Draw();
    }
    ++ it_active_muonid;
  }
  _graphHitPointPhiZAll->Draw("Psame");

  std::map<int,TGraph*>::iterator it_mapClusteredHitPhiZ = _mapClusteredHitPhiZ.begin();
  while( it_mapClusteredHitPhiZ != _mapClusteredHitPhiZ.end() ){
    it_mapClusteredHitPhiZ->second->Draw("Psame");
    ++it_mapClusteredHitPhiZ;
  }

  return;
}

double TrkFinding::readMCInfo(){
  log("debug") << "*****[START: readMCInfo]*****" << std::endl;
  double maximum_positron_energy = 0.0;
  const vector<const MCParticle*> &mc = TreeManager::instance()->getBranchVec<MCParticle>(_parMCParticleName.c_str());

  int muonid;
  for( int ipar=0; ipar<mc.size(); ipar++ ){
    if( mc.at(ipar)->_pdg!=-13 ) continue;
    muonid = mc.at(ipar)->_muonID;
    for( int idau=0; idau<mc.at(ipar)->_daughters.size(); idau++ ){
      if( mc.at(ipar)->_daughters.at(idau)->_pdg==-11 ){
	log("debug") << "ipar = "        << ipar
		     << " : pdg = "      << mc.at(ipar)->_pdg
		     << ", muonID = "    << mc.at(ipar)->_muonID
		     << ", Ndaughter = " << mc.at(ipar)->_daughters.size()
		     << ", idau = "      << idau
		     << " : pdg = "      << mc.at(ipar)->_daughters.at(idau)->_pdg
		     << ", E(p) = "      << mc.at(ipar)->_daughters.at(idau)->_p.E() << " MeV"
		     << std::endl;
	if( muonid>=0 ){
	  _mapMCMuon    [muonid] = mc.at(ipar);
	  _mapMCPositron[muonid] = mc.at(ipar)->_daughters.at(idau);
	  std::vector<int> tmp_vector;
	  _mapMCTrack[muonid] = tmp_vector;
	}
	if( maximum_positron_energy < mc.at(ipar)->_daughters.at(idau)->_p.E() ) maximum_positron_energy = mc.at(ipar)->_daughters.at(idau)->_p.E();
	break;
      }
    }
  }
  log("debug") << " ==> " << _mapMCPositron.size() << " muons are found in total." << std::endl;

  int index;
  for( int ivec=0; ivec<_hitInfoRecoHits.size(); ivec++ ){
    index = _orderMCTruthTime.at(ivec);
    muonid = _hitInfoMuonID.at(index);
    if( _hitInfoRecoHits.at(index)->_trueType ) continue; // remove ghost or noise hit.
    if( muonid>=0 ){
      _mapMCTrack[muonid].push_back(index);
    }
  }

  return maximum_positron_energy;
}

void TrkFinding::drawClusteringResult(){
  if( !_parDrawLevel ) return;
  log("debug") << "*****[START : drawClusteringResult]*****" << std::endl;
  log("debug") << "current maximum clusterNo is " << _clusterNoMax << std::endl;
  int already_plotted_clusterNo = -999;

  if( _mapClusteredHitXY.size() ){
    std::map<int,TGraph*>::iterator it = _mapClusteredHitXY.end();
    --it;
    log("debug") << "TGraph for clustered hit has been already plotted up to clusterNo = " << it->first << std::endl;
    already_plotted_clusterNo = it->first;
  }else{
    log("debug") << "TGraph for clustered hit has not been plotted at all" << std::endl;
  }

  // hough transformation
  _canvas1->cd(3);
  _houghPhiZHist.at(_houghPhiZHist.size()-1)->Draw("COLZ");

  // hough fit
  _canvas1->cd(5);
  std::map<int,TF1*>::iterator it_houghPhiZFunc = _houghPhiZFunc.begin();
  while(it_houghPhiZFunc != _houghPhiZFunc.end()){
    int  lineNo = it_houghPhiZFunc->first;
    TF1* func   = it_houghPhiZFunc->second;
    func->SetLineColor(lineNo%_cycleColor+2);
    func->Draw("Lsame");
    ++it_houghPhiZFunc;
  }

  // clustered hit
  std::map<int,int>::iterator it_cluster_line = _mapActiveClusterNoLineNo.begin();
  int prev_lineNo;
  int cnt_col_cluster=0;
  while(it_cluster_line != _mapActiveClusterNoLineNo.end()){
    int clusterNo = it_cluster_line->first;
    int lineNo    = it_cluster_line->second;
    if( it_cluster_line==_mapActiveClusterNoLineNo.begin() ) prev_lineNo = lineNo;

    if( clusterNo <= already_plotted_clusterNo ){

      if( _mapClusteredHitXY[clusterNo]->GetN()!=_mapClusteringResult[clusterNo].size() ){ // Hits might be added into clustered due to new hits in new time bin
	for( int ip=_mapClusteredHitXY[clusterNo]->GetN(); ip<_mapClusteringResult[clusterNo].size(); ip++ ){
	  int index = _mapClusteringResult[clusterNo][ip];
	  _mapClusteredHitXY  [clusterNo]->SetPoint( _mapClusteredHitXY  [clusterNo]->GetN(),_hitInfoRecoHits.at(index)->_pos.X(),_hitInfoRecoHits.at(index)->_pos.Y()  );
	  _mapClusteredHitPhiZ[clusterNo]->SetPoint( _mapClusteredHitPhiZ[clusterNo]->GetN(),_hitInfoPhi.at(index),               _hitInfoRecoHits.at(index)->_pos.Z()  );
	}
      }


      
      ++it_cluster_line;
      continue;
    }
    
    if( prev_lineNo!=lineNo ) cnt_col_cluster = 0;
    
    TGraph* g_clustered_hit_xy   = makeGraphClusterNoXY  (clusterNo);
    TGraph* g_clustered_hit_phiz = makeGraphClusterNoPhiZ(clusterNo);
    g_clustered_hit_xy  ->SetMarkerColor(_mapActiveClusterNoLineNo[clusterNo]%_cycleColor+2);
    g_clustered_hit_xy  ->SetMarkerStyle(20+cnt_col_cluster);
    g_clustered_hit_xy  ->SetMarkerSize(0.7);
    g_clustered_hit_phiz->SetMarkerColor(_mapActiveClusterNoLineNo[clusterNo]%_cycleColor+2);
    g_clustered_hit_phiz->SetMarkerStyle(20+cnt_col_cluster);
    g_clustered_hit_phiz->SetMarkerSize(0.7);

    _mapClusteredHitXY  [clusterNo] = g_clustered_hit_xy;
    _mapClusteredHitPhiZ[clusterNo] = g_clustered_hit_phiz;
    _canvas1->cd(4);
    g_clustered_hit_xy->Draw("Psame");
    _canvas1->cd(5);
    g_clustered_hit_phiz->Draw("Psame");
    log("debug") << " => draw hit cluster with clusterNo = " << clusterNo << ", lineNo = " << lineNo << std::endl;

    ++it_cluster_line;
    cnt_col_cluster++;
    prev_lineNo = lineNo;
  }

  return;
}

int TrkFinding::getMuonIDfromRecoHit( const RecoHit* recohit ){
  if( recohit->_stripClusters.size()==0 ) return -999;
  if( recohit->_stripClusters.at(0)->_stripHits.size()==0 ) return -999;
  if( recohit->_stripClusters.at(0)->_stripHits.at(0)->_simHits.size()==0 ) return -999;
  if( !recohit->_stripClusters.at(0)->_stripHits.at(0)->_simHits.at(0)->getMCStep() ) return -999;
  if( !recohit->_stripClusters.at(0)->_stripHits.at(0)->_simHits.at(0)->getMCStep()->getMCP() ) return -999;
  return recohit->_stripClusters.at(0)->_stripHits.at(0)->_simHits.at(0)->getMCStep()->getMCP()->_muonID;
}

double TrkFinding::getMCTruthTimefromRecoHit( const RecoHit* recohit ){
  if( recohit->_stripClusters.size()==0 ) return -999;
  if( recohit->_stripClusters.at(0)->_stripHits.size()==0 ) return -999;
  if( recohit->_stripClusters.at(0)->_stripHits.at(0)->_simHits.size()==0 ) return -999;
  if( !recohit->_stripClusters.at(0)->_stripHits.at(0)->_simHits.at(0)->getMCStep() ) return -999;
  return recohit->_stripClusters.at(0)->_stripHits.at(0)->_simHits.at(0)->getMCStep()->_time;
}

void TrkFinding::checkRecoHits( const RecoHit* recohit ){
  log("debug") << "*****[START : checkRecoHits]*****" << std::endl;
  log("debug") << "   TrueType = " << recohit->_trueType << " : N_stripclsuters = " << recohit->_stripClusters.size() << std::endl;
  for( int istripclusters=0; istripclusters<recohit->_stripClusters.size(); istripclusters++ ){
    //log("debug") << "      strip-ID = " << recohit->_stripClusters.at(istripclusters)->_stripID
    //<< ", N_striphits = "  << recohit->_stripClusters.at(istripclusters)->_stripHits.size() << std::endl;
    for( int istriphits=0; istriphits<recohit->_stripClusters.at(istripclusters)->_stripHits.size(); istriphits++ ){
      //log("debug") << "         strip-ID = " << recohit->_stripClusters.at(istripclusters)->_stripHits.at(istriphits)->_stripID
      //<< ", N_simhits = " << recohit->_stripClusters.at(istripclusters)->_stripHits.at(istriphits)->_simHits.size() << std::endl;
      for( int isimhits=0; isimhits<recohit->_stripClusters.at(istripclusters)->_stripHits.at(istriphits)->_simHits.size(); isimhits++ ){
	if( recohit->_stripClusters.at(istripclusters)->_stripHits.at(istriphits)->_simHits.at(isimhits)->getMCStep() ){
	  if( recohit->_stripClusters.at(istripclusters)->_stripHits.at(istriphits)->_simHits.at(isimhits)->getMCStep()->getMCP() ){
	    const MCParticle* mcp = recohit->_stripClusters.at(istripclusters)->_stripHits.at(istriphits)->_simHits.at(isimhits)->getMCStep()->getMCP();
	    //log("debug") << "            pdgCode = " << recohit->_stripClusters.at(istripclusters)->_stripHits.at(istriphits)->_simHits.at(isimhits)->getMCStep()->getMCP()->_pdg
	    //<< ", muon-ID = "           << recohit->_stripClusters.at(istripclusters)->_stripHits.at(istriphits)->_simHits.at(isimhits)->getMCStep()->getMCP()->_muonID
	    //<<" , time(mcp) = "         << recohit->_stripClusters.at(istripclusters)->_stripHits.at(istriphits)->_simHits.at(isimhits)->getMCStep()->getMCP()->_time
	    //<<" , time(mcstep) = "      << recohit->_stripClusters.at(istripclusters)->_stripHits.at(istriphits)->_simHits.at(isimhits)->getMCStep()->_time
	    //<<" , time(sim) = "         << recohit->_stripClusters.at(istripclusters)->_stripHits.at(istriphits)->_simHits.at(isimhits)->_time
	    //<< std::endl;
	    ;
	  }
	}
      }
    }
  }
  return;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int TrkFinding::clusteringNextTimeBin( double time_window_min, double time_window_max ){
  log("debug") << "*****[START : clusteringNextTimeBin]*****" << std::endl;
  log("debug") << "Ncluster = "           << getNClusters()
	       << ", Ncluster(active) = " << getNActiveClusters()
	       << std::endl;

  std::map<int,int>::iterator it_active_clusterNo = _mapActiveClusterNo.begin();
  while( it_active_clusterNo != _mapActiveClusterNo.end() ){ // Begin of Loop for active-clusterNo
    int               clusterNo    = it_active_clusterNo->first; // clusterNo
    std::vector<int>& extrap_index = _mapClusteringResult[clusterNo];

    int prev_nhit = extrap_index.size();
    expandClustering( extrap_index, true  ); // expand clustering along only forward direction
    if( prev_nhit == extrap_index.size() ){
      ++it_active_clusterNo;
      continue;
    }

    log("debug") << extrap_index.size() - prev_nhit << " hits are added to clusterNo = " << clusterNo << std::endl;
    int cnt_seq = prev_nhit;

    // Update clustering result if new hits are added.
    for( unsigned int ihit=cnt_seq; ihit<extrap_index.size(); ihit++ ){
      log("debug") << "index = "       << extrap_index.at(ihit)
		   << ", vaneID = "    << _hitInfoVaneID.at(extrap_index.at(ihit))
		   << ", clusterNo = " << clusterNo
		   << std::endl;
      _hitInfoClusterNo.at ( extrap_index.at(ihit) ) = clusterNo;
      _hitInfoSequenceNo.at( extrap_index.at(ihit) ) = cnt_seq;
      cnt_seq++;
    }
    
    _mapActiveClusterNo[clusterNo] = _timeWindowMax;

    ++it_active_clusterNo;    
  } // End of Loop for active-clusterNo
  log("debug") << "   *****[FINISH : Clustering_NextTimeBin]*****" << std::endl;
  return 1;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TrkFinding::makeGraphXY( TGraph* g_all, TGraph* g_sig, TGraph* g_bkg, double time_window_min, double time_window_max ){
  if( time_window_min<0 ) time_window_min = _timeWindowMin;
  if( time_window_max<0 ) time_window_max = _timeWindowMax;
  makeGraph( g_all, g_sig, g_bkg, 0, time_window_min, time_window_max );
  return;
}
void TrkFinding::makeGraphPhiZ( TGraph* g_all, TGraph* g_sig, TGraph* g_bkg, double time_window_min, double time_window_max ){
  if( time_window_min<0 ) time_window_min = _timeWindowMin;
  if( time_window_max<0 ) time_window_max = _timeWindowMax;
  makeGraph( g_all, g_sig, g_bkg, 1, time_window_min, time_window_max );
  return;
}
void TrkFinding::makeGraph( TGraph* g_all, TGraph* g_sig, TGraph* g_bkg, int sel, double time_window_min, double time_window_max ){
  g_all->Set(0); // dangerous ???
  g_sig->Set(0); // dangerous ???
  g_bkg->Set(0); // dangerous ???
  
  if( time_window_min<0 ) time_window_min = _timeWindowMin;
  if( time_window_max<0 ) time_window_max = _timeWindowMax;

  for( int ihit=0; ihit<_hitInfoRecoHits.size(); ihit++ ){
    if( time_window_min > _hitInfoRecoHits.at(ihit)->_time || time_window_max < _hitInfoRecoHits.at(ihit)->_time ) continue; // time-window is "min <= time <= max"
    int muonid = getMuonIDfromRecoHit( _hitInfoRecoHits.at(ihit));
    //int energy = (muonid >=0 ? _mapMCPositron[getMuonIDfromRecoHit( _hitInfoRecoHits.at(ihit))]->_p.E() : 0.0);
    int momentum = (muonid >=0 ? _mapMCPositron[getMuonIDfromRecoHit( _hitInfoRecoHits.at(ihit))]->_p.P() : 0.0);
    if( sel==0 ){
      g_all->SetPoint( g_all->GetN(), _hitInfoRecoHits.at(ihit)->_pos.X(), _hitInfoRecoHits.at(ihit)->_pos.Y() );
      //if( energy > _parEvtDispThresholdSignalEnergy ) g_sig->SetPoint( g_sig->GetN(), _hitInfoRecoHits.at(ihit)->_pos.X(), _hitInfoRecoHits.at(ihit)->_pos.Y() );
      if( momentum > _parEvtDispThresholdMomentum ) g_sig->SetPoint( g_sig->GetN(), _hitInfoRecoHits.at(ihit)->_pos.X(), _hitInfoRecoHits.at(ihit)->_pos.Y() );
      else                                            g_bkg->SetPoint( g_bkg->GetN(), _hitInfoRecoHits.at(ihit)->_pos.X(), _hitInfoRecoHits.at(ihit)->_pos.Y() );
    }else if( sel==1 ){
      g_all->SetPoint( g_all->GetN(), _hitInfoPhi.at(ihit), _hitInfoRecoHits.at(ihit)->_pos.Z() );
      //if( energy > _parEvtDispThresholdSignalEnergy ) g_sig->SetPoint( g_sig->GetN(), _hitInfoPhi.at(ihit),                _hitInfoRecoHits.at(ihit)->_pos.Z() );
      if( momentum > _parEvtDispThresholdMomentum ) g_sig->SetPoint( g_sig->GetN(), _hitInfoPhi.at(ihit),                _hitInfoRecoHits.at(ihit)->_pos.Z() );
      else                                            g_bkg->SetPoint( g_bkg->GetN(), _hitInfoPhi.at(ihit),                _hitInfoRecoHits.at(ihit)->_pos.Z() );
    }else{
      log("debug") << "[ABORT] Wrong selection : " << sel << std::endl, abort();
    }
  }
  
  return;
}

TGraph* TrkFinding::makeGraphClusterNoXY( int clusterNo, double time_window_min, double time_window_max ){
  if( time_window_min<0 ) time_window_min = _timeWindowMin;
  if( time_window_max<0 ) time_window_max = _timeWindowMax;
  return makeGraphClusterNo( 0, clusterNo, time_window_min, time_window_max );
}
TGraph* TrkFinding::makeGraphClusterNoPhiZ( int clusterNo, double time_window_min, double time_window_max ){
  if( time_window_min<0 ) time_window_min = _timeWindowMin;
  if( time_window_max<0 ) time_window_max = _timeWindowMax;
  return makeGraphClusterNo( 1, clusterNo, time_window_min, time_window_max );
}
TGraph* TrkFinding::makeGraphClusterNo( int sel, int clusterNo, double time_window_min, double time_window_max ){
  if( time_window_min<0 ) time_window_min = _timeWindowMin;
  if( time_window_max<0 ) time_window_max = _timeWindowMax;

  if( clusterNo>0 && _mapClusteringResult.count(clusterNo)==0 ) log("debug") << "[ABORT] Wrong clusterNo : " << clusterNo << std::endl, abort();

  TGraph* g = new TGraph();
  if( clusterNo<0 ){ // not-clustered-hits
    for( int ihit=0; ihit<_hitInfoRecoHits.size(); ihit++ ){
      if( time_window_min > _hitInfoRecoHits.at(ihit)->_time || time_window_max < _hitInfoRecoHits.at(ihit)->_time ) continue;
      if( _hitInfoClusterNo.at(ihit) >= 0 ) continue;
      if     ( sel==0 ) g->SetPoint( g->GetN(), _hitInfoRecoHits.at(ihit)->_pos.X(), _hitInfoRecoHits.at(ihit)->_pos.Y() );
      else if( sel==1 ) g->SetPoint( g->GetN(), _hitInfoPhi.at(ihit),                _hitInfoRecoHits.at(ihit)->_pos.Z() );
      else log("debug") << "[ABORT] Wrong selection : " << sel << std::endl, abort();
    }
  }else{ // clustered-hits
    std::vector<int>& index_list = _mapClusteringResult.at(clusterNo);
    for( int ihit=0; ihit<index_list.size(); ihit++ ){
      if     ( sel==0 ) g->SetPoint( g->GetN(), _hitInfoRecoHits.at(index_list.at(ihit))->_pos.X(), _hitInfoRecoHits.at(index_list.at(ihit))->_pos.Y() );
      else if( sel==1 ) g->SetPoint( g->GetN(), _hitInfoPhi.at(index_list.at(ihit)),               _hitInfoRecoHits.at(index_list.at(ihit))->_pos.Z() );
      else log("debug") << "[ABORT] Wrong selection : " << sel << std::endl, abort();
    }
  }
  
  return g;
}

TGraph* TrkFinding::makeGraphLineNoXY( int lineNo, double time_window_min, double time_window_max ){
  if( time_window_min<0 ) time_window_min = _timeWindowMin;
  if( time_window_max<0 ) time_window_max = _timeWindowMin;
  return makeGraphLineNo( 0, lineNo, time_window_min, time_window_max );
}
TGraph* TrkFinding::makeGraphLineNoPhiZ( int lineNo, double time_window_min, double time_window_max ){
  if( time_window_min<0 ) time_window_min = _timeWindowMin;
  if( time_window_max<0 ) time_window_max = _timeWindowMax;
  return makeGraphLineNo( 1, lineNo, time_window_min, time_window_max );
}
TGraph* TrkFinding::makeGraphLineNo( int sel, int lineNo, double time_window_min, double time_window_max ){
  if( time_window_min<0 ) time_window_min = _timeWindowMin;
  if( time_window_max<0 ) time_window_max = _timeWindowMax;

  TGraph* g = new TGraph();
  for( int ihit=0; ihit<_hitInfoRecoHits.size(); ihit++ ){
    if( time_window_min > _hitInfoRecoHits.at(ihit)->_time || time_window_max < _hitInfoRecoHits.at(ihit)->_time ) continue;
    if( _hitInfoLineNoPhiZ.at(ihit) != lineNo                              ) continue;
    if     ( sel==0 ) g->SetPoint( g->GetN(), _hitInfoRecoHits.at(ihit)->_pos.X(), _hitInfoRecoHits.at(ihit)->_pos.Y() );
    else if( sel==1 ) g->SetPoint( g->GetN(), _hitInfoPhi.at(ihit),                _hitInfoRecoHits.at(ihit)->_pos.Z() );
    else log("debug") << "[ABORT] Wrong selection : " << sel << std::endl, abort();
  }
  
  return g;
}

TGraph* TrkFinding::makeGraphMuonIDXY( int muonID, double time_window_min, double time_window_max ){
  if( time_window_min<0 ) time_window_min = _timeWindowMin;
  if( time_window_max<0 ) time_window_max = _timeWindowMin;
  return makeGraphMuonID( 0, muonID, time_window_min, time_window_max );
}
TGraph* TrkFinding::makeGraphMuonIDPhiZ( int muonID, double time_window_min, double time_window_max ){
  if( time_window_min<0 ) time_window_min = _timeWindowMin;
  if( time_window_max<0 ) time_window_max = _timeWindowMax;
  return makeGraphMuonID( 1, muonID, time_window_min, time_window_max );
}
TGraph* TrkFinding::makeGraphMuonID( int sel, int muonID, double time_window_min, double time_window_max ){
  if( time_window_min<0 ) time_window_min = _timeWindowMin;
  if( time_window_max<0 ) time_window_max = _timeWindowMax;

  TGraph* g = new TGraph();
  int index;
  for( int ihit=0; ihit<_hitInfoRecoHits.size(); ihit++ ){
    index = _orderMCTruthTime.at(ihit);
    if( time_window_min > _hitInfoRecoHits.at(index)->_time ||
	time_window_max < _hitInfoRecoHits.at(index)->_time ) continue;
    if( _hitInfoMuonID.at(index) != muonID                  ) continue;
    if( _hitInfoRecoHits.at(index)->_trueType               ) continue; // remove ghost or noise hit.
    if     ( sel==0 ) g->SetPoint( g->GetN(), _hitInfoRecoHits.at(index)->_pos.X(), _hitInfoRecoHits.at(index)->_pos.Y() );
    else if( sel==1 ) g->SetPoint( g->GetN(), _hitInfoPhi.at(index),                _hitInfoRecoHits.at(index)->_pos.Z() );
    else log("debug") << "[ABORT] Wrong selection : " << sel << std::endl, abort();
  }
  
  return g;
}

void TrkFinding::makeMuonDecayPlot(){
  const vector<const MCParticle*> &mc = TreeManager::instance()->getBranchVec<MCParticle>(_parMCParticleName.c_str());
  std::map<int,const MCParticle*>::iterator it_muon     = _mapMCMuon.begin();
  std::map<int,const MCParticle*>::iterator it_positron = _mapMCPositron.begin();
  while(it_muon != _mapMCMuon.end() ||
	it_positron != _mapMCPositron.end()
	){
    const MCParticle* mc_muon     = it_muon->second;
    const MCParticle* mc_positron = it_positron->second;
    int muonid                    = it_muon->first;
    if( muonid != it_positron->first ) log("debug") << "[WARNING] Wrong map : " << muonid << " : " << it_positron->first  << std::endl;

    if( muonid>=0 ){
      double dec_x  = mc_positron->_prodVertex.X();
      double dec_y  = mc_positron->_prodVertex.Y();
      double dec_z  = mc_positron->_prodVertex.Z();
      double dec_px = mc_positron->_p.X();
      double dec_py = mc_positron->_p.Y();
      double dec_pz = mc_positron->_p.Z();
      
      TMarker* decpoint_xy = new TMarker( dec_x, dec_y, 20 );
      decpoint_xy->SetMarkerStyle(20);
      decpoint_xy->SetMarkerSize(0.7);
      _mapDecayPointXY[muonid] = decpoint_xy;
      
      TMarker* decpoint_phiz = new TMarker( calcPhi(dec_y, dec_x), dec_z, 20 );
      decpoint_phiz->SetMarkerColor(2);
      decpoint_phiz->SetMarkerStyle(20);
      decpoint_phiz->SetMarkerSize(0.7);
      _mapDecayPointPhiZ[muonid] = decpoint_phiz; 
      
      TArrow* g_decvec_xy = new TArrow( dec_x, dec_y,
					dec_x+0.5*dec_px, dec_y+0.5*dec_py,
					0.01,">"
					);
      
      TArrow* g_decvec_phiz = new TArrow( calcPhi(dec_y,dec_x),
					  dec_z,
					  calcPhi(dec_y,dec_x) + ( calcPhi(dec_y+0.01*dec_py,dec_x+0.01*dec_px) - calcPhi(dec_y,dec_x) )*100,
					  dec_z + 1.0*dec_pz,
					  0.01,">"
					  );
      
      _mapDecayDirectionXY  [muonid] = g_decvec_xy;
      _mapDecayDirectionPhiZ[muonid] = g_decvec_phiz;
      
      //if( mc_positron->_p.E() > _parEvtDispThresholdSignalEnergy ){
      if( mc_positron->_p.P() > _parEvtDispThresholdMomentum ){
	decpoint_xy  ->SetMarkerColor(2);
	decpoint_phiz->SetMarkerColor(2);
	g_decvec_xy  ->SetLineColor(2);
	g_decvec_phiz->SetLineColor(2);
      }else{
	decpoint_xy  ->SetMarkerColor(1);
	decpoint_phiz->SetMarkerColor(1);
	g_decvec_xy  ->SetLineColor(1);
	g_decvec_phiz->SetLineColor(1);
      }
      
      double ideal_r = sqrt( pow(dec_px,2) + pow(dec_py,2) )/0.9;
      double ideal_slope = TMath::Abs(dec_px/dec_py);
      double ideal_x0 = dec_x;
      double ideal_y0 = dec_y;
      if( dec_py > 0 ) ideal_x0 += ideal_r/sqrt( 1 + pow(ideal_slope,2) );
      else             ideal_x0 -= ideal_r/sqrt( 1 + pow(ideal_slope,2) );
      if( dec_px > 0 ) ideal_y0 -= ideal_r*ideal_slope/sqrt( 1 + pow(ideal_slope,2) );
      else             ideal_y0 += ideal_r*ideal_slope/sqrt( 1 + pow(ideal_slope,2) );
      TArc* g_ideal_trk = new TArc( ideal_x0, ideal_y0, ideal_r );
      g_ideal_trk->SetFillStyle(0);
      g_ideal_trk->SetLineColor(kGray);
      g_ideal_trk->SetLineStyle(3);
      _mapIdealTrajectory[muonid] = g_ideal_trk;
    }

    ++it_muon;
    ++it_positron;
  } 
  return;
}

void TrkFinding::expandClustering( std::vector<int> &extrap_index, bool isforward ){
  if( isforward ) log("debug") << "*****[START : expandClustering (forward)]*****"  << std::endl;
  else            log("debug") << "*****[START : expandClustering (backward)]*****" << std::endl;

  if( extrap_index.size()<3 ) return;

  double x0, y0, r, dphi, extrap_r, extrap_z;
  double extrap_r_same,extrap_z_same,extrap_r_incre,extrap_z_incre,extrap_r_decre,extrap_z_decre;
  double dphi_same_vane, dphi_incre_vane, dphi_decre_vane;
  //bool   match;

  int cnt_miss = 0;
  int index1, index2, index3;
  while( cnt_miss < _parCutExtrapMiss ){
    if( isforward ){
      index1 = extrap_index.at(extrap_index.size()-3);
      index2 = extrap_index.at(extrap_index.size()-2);
      index3 = extrap_index.at(extrap_index.size()-1);
    }else{
      index1 = extrap_index.at(2);
      index2 = extrap_index.at(1);
      index3 = extrap_index.at(0);
    }

    if(_hitInfoVaneID.at(index1)==_hitInfoVaneID.at(index2) || _hitInfoVaneID.at(index2)==_hitInfoVaneID.at(index3)) return;
    
    log("debug") << "      Expansion from "
		<< "vaneID=" << _hitInfoVaneID.at( index1 ) << "(index=" << index1 << ") & "
		<< "vaneID=" << _hitInfoVaneID.at( index2 ) << "(index=" << index2 << ") & "
		<< "vaneID=" << _hitInfoVaneID.at( index3 ) << "(index=" << index3 << ")" << std::endl;

    int result;
    int target_VaneID, incre_VaneID, decre_VaneID;
    int fl_incre;
    int cnt_nocrosspoint     = 0;
    int cnt_change_direction = 0;
    int result_same_vane, result_incre_vane, result_decre_vane;
    // determine next target vane-ID and expected point
    if( cnt_miss==0 ) target_VaneID = _hitInfoVaneID.at(index3);
    while(1){
      log("debug") << "************ cnt_miss = " << cnt_miss
		  << ", cnt_nocrosspoint = " << cnt_nocrosspoint
		  << ", cnt_change_direction = " << cnt_change_direction
		  << std::endl;
      //if(_hitInfoVaneID.at(index1)==_hitInfoVaneID.at(index2) || _hitInfoVaneID.at(index2)==_hitInfoVaneID.at(index3)) return;
      incre_VaneID = incrementVaneID(target_VaneID);
      decre_VaneID = decrementVaneID(target_VaneID);
      if( cnt_miss==0 && cnt_nocrosspoint==0 ){
	//result_same_vane  = extrapolation( index1, index2, index3,                 target_VaneID,  extrap_r, extrap_z, x0, y0, r, dphi_same_vane,  false);
	//result_incre_vane = extrapolation( index1, index2, index3, incrementVaneID(target_VaneID), extrap_r, extrap_z, x0, y0, r, dphi_incre_vane);
	//result_decre_vane = extrapolation( index1, index2, index3, decrementVaneID(target_VaneID), extrap_r, extrap_z, x0, y0, r, dphi_decre_vane);
	result_same_vane  = extrapolation( index1, index2, index3, target_VaneID, extrap_r_same, extrap_z_same, x0, y0, r, dphi_same_vane,  false);
	result_incre_vane = extrapolation( index1, index2, index3, incre_VaneID, extrap_r_incre, extrap_z_incre, x0, y0, r, dphi_incre_vane);
	result_decre_vane = extrapolation( index1, index2, index3, decre_VaneID, extrap_r_decre, extrap_z_decre, x0, y0, r, dphi_decre_vane);
	log("debug") << "result(same ) = " << result_same_vane  << ", dphi = " << dphi_same_vane  << std::endl;
	log("debug") << "result(incre) = " << result_incre_vane << ", dphi = " << dphi_incre_vane << std::endl;
	log("debug") << "result(decre) = " << result_decre_vane << ", dphi = " << dphi_decre_vane << std::endl;
	if( TMath::Abs(dphi_same_vane ) < TMath::Abs(dphi_incre_vane) && TMath::Abs(dphi_same_vane ) < TMath::Abs(dphi_decre_vane) ){
	  //target_VaneID = target_VaneID;
	  //result = extrapolation( index1, index2, index3, target_VaneID, extrap_r, extrap_z, x0, y0, r, dphi, false);
	  fl_incre      = 0;
	  extrap_r = extrap_r_same;
	  extrap_z = extrap_z_same;
	  dphi = dphi_same_vane;
	  result = result_same_vane;
	  cnt_change_direction++;
	  log("debug") << " => SAME is selected : Target vane-ID = " << target_VaneID << std::endl;
	  if( !_parDoSmallCurlFinding ) return;
	}else if( TMath::Abs(dphi_incre_vane) < TMath::Abs(dphi_same_vane) && TMath::Abs(dphi_incre_vane) < TMath::Abs(dphi_decre_vane) ){
	  //target_VaneID = incrementVaneID(target_VaneID);
	  //result        = extrapolation( index1, index2, index3, target_VaneID, extrap_r, extrap_z, x0, y0, r, dphi);
	  target_VaneID = incre_VaneID;
	  extrap_r = extrap_r_incre;
	  extrap_z = extrap_z_incre;
	  dphi = dphi_incre_vane;
	  result = result_incre_vane;
	  fl_incre      = 1;
	  log("debug") << " => INCRE is selected : Target vane-ID = " << target_VaneID << std::endl;
	}else if( TMath::Abs(dphi_decre_vane) < TMath::Abs(dphi_same_vane ) && TMath::Abs(dphi_decre_vane) < TMath::Abs(dphi_incre_vane ) ){
	  //target_VaneID = decrementVaneID(target_VaneID);
	  //result        = extrapolation( index1, index2, index3, target_VaneID, extrap_r, extrap_z, x0, y0, r, dphi);
	  target_VaneID = decre_VaneID;
	  extrap_r = extrap_r_decre;
	  extrap_z = extrap_z_decre;
	  dphi = dphi_decre_vane;
	  result = result_decre_vane;
	  fl_incre      = -1;
	  log("debug") << " => DECRE is selected : Target vane-ID = " << target_VaneID << std::endl;
	}else{
	  log("debug") << " => return due to all result < 0" << std::endl;
	  return;	  
	}
	extrapolation( index1, index2, index3, target_VaneID, extrap_r, extrap_z, x0, y0, r, dphi );
      }else{ // cnt_miss!=0 || cnt_nocrosspoint!=0
	if( fl_incre==0 ){ // for small curl track
	  /*
	  log("debug") << "isMoveAwayFromOrigin(same:incre) : " << isMoveAwayFromOrigin( 1,incrementVaneID(target_VaneID), _hitInfoVaneID.at(index3)) << std::endl;
	  log("debug") << "isMoveAwayFromOrigin(same:decre) : " << isMoveAwayFromOrigin(-1,decrementVaneID(target_VaneID), _hitInfoVaneID.at(index3)) << std::endl;
	  result_incre_vane = extrapolation( index1, index2, index3, incrementVaneID(target_VaneID), extrap_r, extrap_z, x0, y0, r, dphi_incre_vane, isMoveAwayFromOrigin( 1,incrementVaneID(target_VaneID), _hitInfoVaneID.at(index3)) );
	  result_decre_vane = extrapolation( index1, index2, index3, decrementVaneID(target_VaneID), extrap_r, extrap_z, x0, y0, r, dphi_decre_vane, isMoveAwayFromOrigin(-1,decrementVaneID(target_VaneID), _hitInfoVaneID.at(index3)) );
	  */
	  bool isMoveAwayFromOrigin_incre = isMoveAwayFromOrigin( 1,incre_VaneID, _hitInfoVaneID.at(index3));
	  bool isMoveAwayFromOrigin_decre = isMoveAwayFromOrigin(-1,decre_VaneID, _hitInfoVaneID.at(index3));
	  log("debug") << "isMoveAwayFromOrigin(same:incre) : " << isMoveAwayFromOrigin_incre << std::endl;
	  log("debug") << "isMoveAwayFromOrigin(same:decre) : " << isMoveAwayFromOrigin_decre << std::endl;
	  result_incre_vane = extrapolation( index1, index2, index3, incre_VaneID, extrap_r_incre, extrap_z_incre, x0, y0, r, dphi_incre_vane, isMoveAwayFromOrigin_incre );
	  result_decre_vane = extrapolation( index1, index2, index3, decre_VaneID, extrap_r_decre, extrap_z_decre, x0, y0, r, dphi_decre_vane, isMoveAwayFromOrigin_decre );
	  if( result_incre_vane > 0 ){
	    fl_incre      =  1;
	    //target_VaneID = incrementVaneID(target_VaneID);
	    //result        = extrapolation( index1, index2, index3, target_VaneID, extrap_r, extrap_z, x0, y0, r, dphi_incre_vane, isMoveAwayFromOrigin( 1,target_VaneID, _hitInfoVaneID.at(index3)) );
	    target_VaneID = incre_VaneID;
	    extrap_r = extrap_r_incre;
	    extrap_z = extrap_z_incre;
	    dphi = dphi_incre_vane;
	    result = result_incre_vane;
	    log("debug") << "change to incre from same: Target vane-ID = " << target_VaneID << std::endl;
	  }else if( result_decre_vane > 0 ){
	    fl_incre      = -1;
	    //target_VaneID = decrementVaneID(target_VaneID);
	    //result        = extrapolation( index1, index2, index3, target_VaneID, extrap_r, extrap_z, x0, y0, r, dphi_decre_vane, isMoveAwayFromOrigin(-1,target_VaneID, _hitInfoVaneID.at(index3)) );
	    target_VaneID = decre_VaneID;
	    extrap_r = extrap_r_decre;
	    extrap_z = extrap_z_decre;
	    dphi = dphi_decre_vane;
	    result = result_decre_vane;
	    log("debug") << "change to decre from same: Target vane-ID = " << target_VaneID << std::endl;
	  }else{
	    log("debug") << " => return due to out of z" << std::endl;
	    return;
	  }
	}else if( fl_incre==1 ){ // increment vane-ID
	  bool isMoveAwayFromOrigin_incre = isMoveAwayFromOrigin(fl_incre, incre_VaneID, _hitInfoVaneID.at(index3));
	  //log("debug") << "isMoveAwayFromOrigin(incre) : " << isMoveAwayFromOrigin(fl_incre,incrementVaneID(target_VaneID), _hitInfoVaneID.at(index3)) << std::endl;
	  //result = extrapolation( index1, index2, index3, incrementVaneID(target_VaneID), extrap_r, extrap_z, x0, y0, r, dphi, isMoveAwayFromOrigin(fl_incre,incrementVaneID(target_VaneID), _hitInfoVaneID.at(index3)) );
	  log("debug") << "isMoveAwayFromOrigin(incre) : " << isMoveAwayFromOrigin_incre << std::endl;
	  result = extrapolation( index1, index2, index3, incre_VaneID, extrap_r, extrap_z, x0, y0, r, dphi, isMoveAwayFromOrigin_incre );
	  if( result > 0 ){
	    //target_VaneID = incrementVaneID(target_VaneID);
	    target_VaneID = incre_VaneID;
	    log("debug") << "keep incre : Target vane-ID = " << target_VaneID << std::endl;
	  }else if( result == -2 ){
	    return;
	  }else{
	    fl_incre = 0;
	    result   = extrapolation( index1, index2, index3, target_VaneID, extrap_r, extrap_z, x0, y0, r, dphi, false);
	    log("debug") << "change to same from incre: Target vane-ID = " << target_VaneID << std::endl;
	    cnt_change_direction++;
	    if( !_parDoSmallCurlFinding ) return;
	  }
	}else if( fl_incre==-1 ){ // decrement vane-ID
	  bool isMoveAwayFromOrigin_decre = isMoveAwayFromOrigin(fl_incre,decrementVaneID(target_VaneID), _hitInfoVaneID.at(index3));
	  //log("debug") << "isMoveAwayFromOrigin(decre) : " << isMoveAwayFromOrigin(fl_incre,decrementVaneID(target_VaneID), _hitInfoVaneID.at(index3)) << std::endl;
	  //result = extrapolation( index1, index2, index3, decrementVaneID(target_VaneID), extrap_r, extrap_z, x0, y0, r, dphi, isMoveAwayFromOrigin(fl_incre,decrementVaneID(target_VaneID), _hitInfoVaneID.at(index3)) );
	  log("debug") << "isMoveAwayFromOrigin(decre) : " << isMoveAwayFromOrigin_decre << std::endl;
	  result = extrapolation( index1, index2, index3, decre_VaneID, extrap_r, extrap_z, x0, y0, r, dphi, isMoveAwayFromOrigin_decre );
	  if( result > 0 ){
	    //target_VaneID = decrementVaneID(target_VaneID);
	    target_VaneID = decre_VaneID;
	    log("debug") << "keep decre : Target vane-ID = " << target_VaneID << std::endl;
	  }else if( result == -2 ){
	    log("debug") << " => return due to out of z" << std::endl;
	    return;
	  }else{
	    fl_incre = 0;
	    result = extrapolation( index1, index2, index3, target_VaneID,  extrap_r, extrap_z, x0, y0, r, dphi, false);
	    log("debug") << "change to same from decre : Target vane-ID = " << target_VaneID << std::endl;
	    cnt_change_direction++;
	    if( !_parDoSmallCurlFinding ) return;
	  }
	}

      }
      log("debug") << "RESULT = " << result << std::endl;
      log("debug") << "extrap_r = " << extrap_r << ", extrap_z = " << extrap_z << std::endl;

      if( result==1 ){ // cross point exists
	break;
      }else if( result==2 ){ // no cross point
	cnt_nocrosspoint++;
	if( cnt_nocrosspoint>=_parCutExtrapNoCross ){
	  target_VaneID = -999;
	  log("debug") << " => no cross point => return" << std::endl;
	  return;
	}
	log("debug") << " => no cross point => continue" << std::endl;
      }
    } // while(1)
    log("debug") << "Next target vane-ID is " << target_VaneID << std::endl;

    // compare actual hits with expected hit point.
    double dev_min   = 10000;
    int    index_min = -999;
    //for( unsigned int ivane=0; ivane<_orderVaneID.size(); ivane++ ){
    //  int index = _orderVaneID.at(ivane);
    auto range = _multimapVaneID.equal_range(target_VaneID);
    for(auto it=range.first; it!=range.second; it++){
      int index = (*it).second;

      //if( _hitInfoVaneID.at(index)    != target_VaneID ) continue;
      //if( _hitInfoVaneID.at(index)    < target_VaneID ) continue;
      //if( _hitInfoVaneID.at(index)    > target_VaneID ) break;
      if( _hitInfoClusterNo.at(index) >= 0             ) continue;
      if( _timeWindowMin > _hitInfoRecoHits.at(index)->_time || _timeWindowMax < _hitInfoRecoHits.at(index)->_time ) continue; // time-window is "min <= time <= max"

      double dev_r = TMath::Abs( extrap_r - _hitInfoRecoHits.at(index)->_pos.Perp() );
      double dev_z = TMath::Abs( extrap_z - _hitInfoRecoHits.at(index)->_pos.Z()    );
      double dev = sqrt(dev_r*dev_r+dev_z*dev_z);
      if( dev > _parCutExtrapToleranceCoeffDphi*TMath::Abs(dphi) && dev > _parCutExtrapTolerance ) continue;

      if( std::find(extrap_index.begin(), extrap_index.end(), index)!=extrap_index.end() ) continue;

      /*
      bool match = false;
      for( unsigned int l=0; l<extrap_index.size(); l++ ){
	if( index==extrap_index.at(l) ){
	  match = true;
	  break;
	}
      }
      if( match ) continue;
      */
      log("debug") << "            index = " << index << " dev(r) = " << dev_r << ", dev(z) = " << dev_z << ", dev = " << dev << std::endl;
      if( dev_min > dev ){
	dev_min   = dev;
	index_min = index;
      }
    }

    if( dev_min < _parCutExtrapToleranceCoeffDphi*TMath::Abs(dphi) || dev_min < _parCutExtrapTolerance ){
      if( isforward ) extrap_index.push_back(index_min);
      else            extrap_index.insert(extrap_index.begin(), index_min);
      cnt_miss = 0;
      log("debug") << "         => added the hit into the cluster : " << index_min << std::endl;
    }else{
      if( !(extrap_r < _geomRInner || extrap_r > _geomROuter) ) cnt_miss++;
      log("debug") << "         => can not find hit points by extrapolation : dev_min = " << dev_min << ", cnt_miss = " << cnt_miss << std::endl;
    }
  }

  return;
}

int TrkFinding::incrementVaneID( const int vaneID, const int dev ){
  //vaneID += dev;
  //if( vaneID >= _parGeomNVane ) vaneID -= _parGeomNVane;
  //return vaneID;
  int incre_vaneID = vaneID + dev;
  if( incre_vaneID >= _parGeomNVane ) incre_vaneID -= _parGeomNVane;
  return incre_vaneID;
}

int TrkFinding::decrementVaneID( const int vaneID, const int dev ){
  //vaneID -= dev;
  //if( vaneID < 0 ) vaneID += _parGeomNVane;
  //return vaneID;
  int decre_vaneID = vaneID - dev;
  if( decre_vaneID < 0 ) decre_vaneID += _parGeomNVane;
  return decre_vaneID;
}

double TrkFinding::calcDeltaPhi( double phi_start, double phi_end, bool isclockwise ){
  if( isclockwise ){
    while( phi_end > phi_start ) phi_end -= 2.0*TMath::Pi();
  }else{
    while( phi_end < phi_start ) phi_end += 2.0*TMath::Pi();
  }
  double dphi = phi_end - phi_start; // dphi > 0 @anti-clockwise, dphi < 0 @clockwise
  
  return dphi;
}

bool TrkFinding::isMoveAwayFromOrigin( int fl_incre, int target_vaneID, int origin_vaneID ){
  if( target_vaneID==origin_vaneID ){
    return false;
  }else if( fl_incre==1 || fl_incre==-1 ){
    int delta;
    if( fl_incre==1 ) delta = target_vaneID - origin_vaneID;
    else              delta = origin_vaneID - target_vaneID;
    while(delta<0) delta += _parGeomNVane;
    if( delta < _parGeomNVane/2 ) return true;
    else                          return false;
  }else{
    return false; // ????
  }
}

void TrkFinding::calcActiveMuonID( double time_window_min, double time_window_max ){
  if( time_window_min < 0 ) time_window_min = _timeWindowMin;
  if( time_window_max < 0 ) time_window_max = _timeWindowMax;
  _setActiveMuonID.clear();

  int index;
  for( int ivec=0; ivec<_hitInfoRecoHits.size(); ivec++ ){
    index = _orderRecTime.at(ivec);
    _setActiveMuonID.insert(_hitInfoMuonID.at(index));
  }
  return;
}

int TrkFinding::getDominantMuonID( std::multiset<int> & cnt_muonid ){
  /*
  std::multiset<int>::reverse_iterator it_rev_cnt_muonid = cnt_muonid.rbegin();
  int muonid_max = *it_rev_cnt_muonid;
  int muonid     = -999;
  int cnt_max    = 0;
  for( unsigned imuon=0; imuon<=muonid_max; imuon++ ){
    if( cnt_muonid.count(imuon)>cnt_max ){
      muonid  = imuon;
      cnt_max = cnt_muonid.count(imuon);
    }
  }
  */
  std::unordered_map<int,int> muonid_map;
  for(auto it_cnt_muonid : cnt_muonid ){
    if( it_cnt_muonid==-999 ) continue;

    if( muonid_map.find(it_cnt_muonid)==muonid_map.end() ){
      muonid_map[it_cnt_muonid] = 1;
    }else{
      muonid_map[it_cnt_muonid]++;
    }
  }

  int muonid = -999;
  int cnt_max = 0;
  for(auto it_muonid_map : muonid_map ){
    int count = it_muonid_map.second;
    if( count>cnt_max ){
      cnt_max = count;
      muonid = it_muonid_map.first;
    }
  }

  return muonid;
}

void TrkFinding::makeRecTimeHist(double time_window_start, double time_window_end)
{
  int nTimeWindow = ceil((time_window_end - time_window_start)/_parTimeWindowStep);
  _histRecTime = new TH1I("histRecTime", "histRecoTime;Time [ns];Number of RecoHits", nTimeWindow, time_window_start, time_window_start + _parTimeWindowStep*nTimeWindow);
  for(unsigned int iorder=0; iorder<_orderRecTime.size(); iorder++){
    int index = _orderRecTime.at(iorder);
    double time = _hitInfoRecoHits.at(index)->_time;
    _histRecTime->Fill( time );
  }
}

bool TrkFinding::skipTimeWindow()
{
  int binx = _histRecTime->FindBin(_timeWindowMin);
  bool skip = true;
  for(unsigned int ibin=0; ibin<ceil(_parTimeWindowWidth/_parTimeWindowStep); ibin++){
    if( _histRecTime->GetBinContent(binx+ibin) ) skip = false;
  }
  return skip;
}
