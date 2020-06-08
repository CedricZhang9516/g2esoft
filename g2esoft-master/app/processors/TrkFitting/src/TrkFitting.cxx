// TrkFitting.cxx

#include "TreeProc/TreeManager.h"
#include "TreeProc/LogStream.h"
//#include "g2track/g2track.h"
#include "TrkFitting.h"
#include "MyBField.h"
#include "MyFinitePlane.h"

#include <TGeoManager.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TROOT.h>

#include <FieldManager.h>
#include <ConstField.h>
#include <MaterialEffects.h>
#include <TGeoMaterialInterface.h>
#include <KalmanFitter.h>
#include <KalmanFitterRefTrack.h>
#include <KalmanFitterInfo.h>
#include <KalmanFittedStateOnPlane.h>
#include <DAF.h>
#include <EventDisplay.h>
#include <AbsTrackRep.h>
#include <RKTrackRep.h>
#include <MeasuredStateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>
#include <PlanarMeasurement.h>
#include <DetPlane.h>

#include <iostream>
#include <stdexcept>
#include <fstream>

using namespace g2esoft;
using namespace TreeProc;
using namespace std;

static Factory<Processor, TrkFitting> _f("TrkFitting");

TrkFitting::TrkFitting(const char *procName, const char *instName) 
  : Processor(procName, instName),
    _fitter(0),
    _matEff(0),
    _display(0)
{}

TrkFitting::~TrkFitting(){}

void TrkFitting::config(Parameter *param){
  _parInputTrackName    = param->get("InputTrackName",    string("Tracks"));
  _parFittedTrackName = param->get("FittedTrackName", string("FittedTracks")); 
  TreeManager::instance()->registerBranch(_parFittedTrackName.c_str(), &_fittedTracks);
  _parNMaxIter = param->get("NMaxIter", (int)8);
  _parNMinIter = param->get("NMinIter", (int)6);
  _parNMaxFailed = param->get("NMaxFailed", (int)-1);
  _parDRelChi2 = param->get("DRelChi2", (double)0.2);
  _parDPVal    = param->get("DPVal",    (double)1.e-3);
  _parFitterName = param->get("FitterName", string("DafRef"));
  _parDetectorGeometry = param->get("DetectorGeometry", string("Detector.gdml"));
  _parShowEventDisplay = param->get("ShowEventDisplay", (bool)false);
  _parUseStripCluster = param->get("UseStripCluster", (bool)true);
  _parRemoveGhostTrack = param->get("RemoveGhostTrack", (bool)true);
  _parConnectTrack = param->get("ConnectTrack", (bool)true);
  _parMomentumMin = param->get("MomentumMin", (double)150.); // MeV/c
  _parChi2NDFMax = param->get("Chi2NDFMax", (double)100.);
  _parConnectChi2NDFMax = param->get("ConnectChi2NDFMax", (double)3.);
  _parConnectTimeDiff = param->get("ConnectTimeDiff", (double)5.); // ns
  _parConnectSearchTime = param->get("ConnectSearchTime", (double)10.); // ns
  _parGhostTimeDiff = param->get("GhostTimeDiff", (double)2.5);
  _parNVane = param->get("NVane",(int)40);
  _parSensorOriginX = param->get("SensorOriginX",vector<double>(91.975, 191.245));
  _parSensorOriginZ = param->get("SensorOriginZ",vector<double>(0.5,99.77));
  _parPosFirstStrip = param->get("PosFirstStrip",(double)0.840);
  _parStripPitch = param->get("StripPitch",(double)0.190);
  _parSensorDistance = param->get("SensorDistance", (double)3.32);
  _parSensorSize = param->get("SensorSize", (double)98.77);
  _parNStripPerRow = param->get("NStripPerRow", (int)512);
  _parTimeDiv = param->get("TimeDiv", (int)5);
 
  _parMagnetData = param->get("MagnetData", string(""));
  _parUniformField = param->get("UniformField", (bool)true);
  _parInterpolateField = param->get("InterpolateField", (bool)true);

  // init geometry and mag. field
  new TGeoManager("DetectorGeometry", "DetectorGeometry");
  ifstream ifs(_parDetectorGeometry.c_str());
  if(!ifs.is_open()){
    throw std::runtime_error("Detector geometry file: "+ _parDetectorGeometry+" cannot be opend.");
  }
  ifs.close();
  if( TGeoManager::Import(_parDetectorGeometry.c_str())==0 ){
    throw std::runtime_error("Failed to import detector geometry: "+_parDetectorGeometry);
  }
  
  if( _parUniformField ){
    const double BZ = 30; // kGauss
    genfit::FieldManager::getInstance()->init(new genfit::ConstField(0., 0., BZ));
  }else{
    genfit::MyBField *bfield = new genfit::MyBField(_parMagnetData);
    if(_parInterpolateField) bfield->setDoInterpolate(_parInterpolateField);
    genfit::FieldManager::getInstance()->init(bfield);
  }
  
  _matEff = genfit::MaterialEffects::getInstance();
  _matEff->init(new genfit::TGeoMaterialInterface());
  
  if(_parShowEventDisplay){
    gROOT->SetBatch(kFALSE);
    _display = genfit::EventDisplay::getInstance();
  }

  // Genfit options
  const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedAverage;

  if( _parFitterName=="SimpleKalman" ){
    _fitter = new genfit::KalmanFitter(_parNMaxIter, _parDPVal);
    _fitter->setMultipleMeasurementHandling(mmHandling);
  }else if( _parFitterName=="RefKalman" ){
    _fitter = new genfit::KalmanFitterRefTrack(_parNMaxIter, _parDPVal);
    _fitter->setMultipleMeasurementHandling(mmHandling);
  }else if( _parFitterName=="DafSimple" ){
    _fitter = new genfit::DAF(false);
  }else if( _parFitterName=="DafRef" ){
    _fitter = new genfit::DAF();
  }else{
    log("error") << "Undefined FitterName: " << _parFitterName << endl;
    throw std::runtime_error("Undefined FitterName");
  }

  _fitter->setMinIterations(_parNMinIter);
  _fitter->setMaxIterations(_parNMaxIter);
  _fitter->setRelChi2Change(_parDRelChi2);
  _fitter->setMaxFailedHits(_parNMaxFailed);
  //_fitter->setDebugLvl(1);

}


void TrkFitting::process(int nEvent, Parameter *param){
  log("warning") << "TrkFitting process starts for event " << nEvent << endl;

  // clear fitted tracks
  for(unsigned int j=0; j<_fittedTracks.size(); j++){
    delete _fittedTracks[j];
  }
  _fittedTracks.clear();
  _tempTracks.clear();

  const int pdgWanted = -11;
  const double resolution = _parStripPitch*0.1/sqrt(12); // cm
  const double timereso = _parTimeDiv/sqrt(12.); // ns

  // input: tracks
  const vector<const Track *> &tracks = TreeManager::instance()->getBranchVec<Track>(_parInputTrackName.c_str());
  log("debug") << tracks.size() << " Tracks are in this event" << endl;
  const Track*   currentTrack = 0;
  const RecoHit* currentRecoHit = 0;
  const StripCluster* currentStripCluster = 0;
  vector<double> thit;
  vector<double> thitTrue;
  vector<double> trackLastTime;
  vector<genfit::Track *> fitTrackVec;
  
  TF1 *f1 = new TF1("f1", "[0]+[1]*x", 0., 1.);
  f1->FixParameter(1, 1.);
  f1->SetLineColor(kBlack);
  TGraphErrors *grHitTime = new TGraphErrors;
  grHitTime->SetNameTitle("hitTime", "hitTime;Time of flight [ns];Reconstructed hit time [ns]");
  TGraphErrors *grTrueHitTime = new TGraphErrors;
  grTrueHitTime->SetNameTitle("trueHitTime", "trueHitTime;Time of flight [ns];True hit time [ns]");

  for(unsigned int j=0; j<tracks.size(); j++){
    currentTrack = tracks.at(j);
    if( currentTrack==0 ){
      log("debug") << j << " th track is empty" << endl;
      continue;
    }

    unsigned int nHits = currentTrack->_recoHits.size();
    log("debug") << "Number of RecoHits is " << nHits << endl;
    if( nHits<3 ) continue;

    if( currentTrack->_p.Mag()<_parMomentumMin ){
      log("info") << "Initail momentum is  " << currentTrack->_p.Mag() << " MeV/c. Skip TrkFitting for this track" << endl;
      continue;
    }

    currentRecoHit = currentTrack->_recoHits.at(0);
    if( currentRecoHit==0 ){
      continue;
    }

    thit.clear();
    thitTrue.clear();
    grHitTime->Set(0);
    grTrueHitTime->Set(0);

    // track for g2track output
    Track *fittedTrack = new Track();
    fittedTrack->_charge = currentTrack->_charge;
    fittedTrack->_time = currentTrack->_time;

    const MCParticle *mcp = currentTrack->getMCP();
    if(mcp){
      log("debug") << "true particle pdgID = " << mcp->_pdg << ", p = " << mcp->_p.P() << " MeV/c (" << mcp->_p.X() << ", " << mcp->_p.Y() << ", " << mcp->_p.Z() << ")" << endl;
      fittedTrack->_mcp = const_cast<MCParticle *>(mcp);
    }

    // set the initial state (for Kalman filter step 0) Units: cm, GeV/c
    TVector3 pos = currentRecoHit->_pos;
    pos *= 0.1; // mm -> cm
    TVector3 mom = currentTrack->_p;
    mom *= 0.001; // MeV/c -> GeV/c
    const int pdg = pdgWanted;
    
    // get the charge
    //const double charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/3.;

    // approximate covariance (Kalman filter step 1)
    TMatrixDSym covM(6);
    for(int k=0; k<3; ++k){
      covM(k,k) = resolution*resolution;
    }
    for(int k=3; k<6; ++k){
      covM(k,k) = pow(resolution/(sqrt(3)*nHits),2);
    }

    genfit::AbsTrackRep *rep = new genfit::RKTrackRep(pdg);
    genfit::MeasuredStateOnPlane stateMeas(rep);
    stateMeas.setPosMomCov(pos, mom, covM);
    
    // create track
    TVectorD seedState(6);
    TMatrixDSym seedCov(6);
    stateMeas.get6DStateCov(seedState, seedCov);
    genfit::Track *fitTrack = new genfit::Track(rep, seedState, seedCov);

    TVector3 truep;
    TVector3 truepos;
    double timeMin = 1.e6;
    const MCStep *firstMCStep = 0;

    for(unsigned int iHit=0; iHit<nHits; iHit++){
      currentRecoHit = currentTrack->_recoHits.at(iHit);
      if( currentRecoHit==0 ){
	continue;
      }
      fittedTrack->_recoHits.push_back(currentRecoHit);
      
      log("debug") << "RecoHit: x = " << currentRecoHit->_pos.X() << " mm, y = " << currentRecoHit->_pos.Y() << " mm, z = " << currentRecoHit->_pos.Z()  << " mm" << endl;

      if( !_parUseStripCluster ){
	thit.push_back( currentRecoHit->_time );
	setRecoHitPosition(fitTrack, currentRecoHit);
      }
      for(unsigned int iClus=0; iClus<currentRecoHit->_stripClusters.size(); iClus++){
	currentStripCluster = currentRecoHit->_stripClusters.at(iClus);
	if( _parUseStripCluster ){
	  thit.push_back( currentStripCluster->_time );
	  setStripClusterPosition(fitTrack, currentStripCluster);
	}
	 
	double timeTrue = 1.e9;
	for(unsigned int iStrip=0; iStrip<currentStripCluster->_stripHits.size(); iStrip++){	      
	  for(unsigned int iSim=0; iSim<currentStripCluster->_stripHits.at(iStrip)->_simHits.size(); iSim++){
	    const SimHit *simhit = currentStripCluster->_stripHits.at(iStrip)->_simHits.at(iSim);
	    if(simhit){
	      if( simhit->_time<timeTrue ){
		timeTrue = simhit->_time;
	      }
	      
	      // get initial hit true momentum
	      if( iHit==0 ){
	      	const MCStep *mcstep = simhit->getMCStep();
		if(mcstep){
		  if( mcstep->getMCP() ){
		    if(mcstep->getMCP()->_pdg==-11){
		      if(mcstep->_time<timeMin){
			timeMin = mcstep->_time;
			truep = mcstep->_p.Vect();
			truepos = mcstep->_pos;
			firstMCStep = mcstep;
		      }
		    }
		  }
		} // mcstep
	      } // iHit==0 && iClus==0
	    } // simhit
	  } // iSim
	} // iStrip
	thitTrue.push_back( timeTrue );
      } // iClus
    } // iHit

    if(firstMCStep){
      log("info") << "true p = " << truep.Mag() << " MeV/c, px = " << truep.X() << " py = " << truep.Y() << " pz = " << truep.Z() << endl;
      log("info") << "true x = " << truepos.X() << " mm, y = " << truepos.Y() << " mm, z = " << truepos.Z() << " mm" << endl;
      fittedTrack->_mcStep = const_cast<MCStep *>(firstMCStep);
    }
    
    // check
    fitTrack->checkConsistency();

    // do the fit
    _fitter->processTrack(fitTrack);

    // check
    fitTrack->checkConsistency();

    if(_parShowEventDisplay) _display->addEvent(fitTrack);

    // check if fit was succcessful
    genfit::FitStatus *fitstat = fitTrack->getFitStatus(rep);
    bool fitconv = fitstat->isFitConverged();
    double chi2 = fitstat->getChi2();
    double NDF = fitstat->getNdf();
    fittedTrack->_chi2ndf = chi2/NDF;
    
    if(!fitconv){
      log("info") << "Track could not be fitted successfully" << endl;
      delete fittedTrack;
      delete fitTrack;
      continue;
    }

    genfit::TrackPoint* tp = fitTrack->getPointWithMeasurementAndFitterInfo(0, rep);
    if(tp==NULL){
      log("info") << "Track has no TrackPoint with fitterInfo" << endl;
      delete fittedTrack;
      delete fitTrack;
      continue;
    }

    if( _parChi2NDFMax>0 && chi2/NDF > _parChi2NDFMax ){
      log("info") << "Chi2/NDF " << chi2/NDF << " is larger than " << _parChi2NDFMax << endl;
      delete fittedTrack;
      delete fitTrack;
      continue;
    }

    genfit::KalmanFittedStateOnPlane kfsop(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()));
    TVector3 posResult = kfsop.getPos();
    posResult *= 10.; // cm -> mm
    TVector3 momResult = kfsop.getMom();
    momResult *= 1000.; // GeV/c -> MeV/c
    TMatrixDSym covResult = kfsop.getCov(); // q/p, u', v', u, v
    log("info") << "fitted p = " << momResult.Mag() << " MeV/c, px = " << momResult.X() << " py = " << momResult.Y() << " pz = " << momResult.Z() << endl;
    log("info") << "fitted x = " << posResult.X() << " mm, y = " << posResult.Y() << " mm, z = " << posResult.Z() << " mm, t = " << kfsop.getTime() << " ns" << endl;
    fittedTrack->_p = momResult;
    fittedTrack->_pos = posResult;
    fittedTrack->_covM.ResizeTo(covResult);
    fittedTrack->_covM = covResult;

    // first hit time calculation
    double tmin=1.e9, tmax=-1.e9;
    int nHitTime=0;
    for(int k=0; k<thit.size(); k++){
      genfit::TrackPoint *tpTemp = fitTrack->getPointWithMeasurementAndFitterInfo(k, rep);
      if(tpTemp==NULL){
	continue;
      }
      genfit::KalmanFittedStateOnPlane kfsopTemp(*(static_cast<genfit::KalmanFitterInfo*>(tpTemp->getFitterInfo(rep))->getBackwardUpdate()));
      //genfit::KalmanFittedStateOnPlane kfsopTemp(*(static_cast<genfit::KalmanFitterInfo*>(tpTemp->getFitterInfo(rep))->getForwardUpdate()));
      double timeTemp = kfsopTemp.getTime();
      if(timeTemp<tmin) tmin = timeTemp;
      if(timeTemp>tmax) tmax = timeTemp;
      grHitTime->SetPoint(nHitTime, timeTemp, thit.at(k) );
      grHitTime->SetPointError(nHitTime, 0, timereso);
      grTrueHitTime->SetPoint(nHitTime, timeTemp, thitTrue.at(k) );
      nHitTime++;
    }

    f1->SetRange(tmin, tmax); 
    grHitTime->Fit("f1","RQ");
    fittedTrack->_time = f1->GetParameter(0);
    log("debug") << "tmin = " << tmin << " tmax = " << tmax << " decaytime = " << f1->GetParameter(0) << endl;

    _tempTracks.push_back(fittedTrack);
    trackLastTime.push_back(tmax+f1->GetParameter(0));
    fitTrackVec.push_back(fitTrack);
  }

  delete grTrueHitTime;
  delete grHitTime;
  delete f1;

  Track* tempTrack;

  //---------------------//
  // ghost track removal //
  //---------------------//
  if( _parRemoveGhostTrack ){
    log("debug") << "[Start ghost track removal]" << endl;
    for(unsigned int itrack=0; itrack<_tempTracks.size(); itrack++){
      tempTrack = _tempTracks.at(itrack);
      double ghostRatio = 0.;
      double stripTimeDiff = 0.;
      for(unsigned int ireco=0; ireco<tempTrack->_recoHits.size(); ireco++){
	currentRecoHit = tempTrack->_recoHits.at(ireco);
	double stripTimeMin = 1e9;
	double stripTimeMax = -1e9;
	for(unsigned int iclus=0; iclus<currentRecoHit->_stripClusters.size(); iclus++){
	  const StripCluster *stripClus = currentRecoHit->_stripClusters.at(iclus);
	  double stripTime = stripClus->_time;
	  if( stripTime<stripTimeMin ) stripTimeMin = stripTime;
	  if( stripTime>stripTimeMax ) stripTimeMax = stripTime;
	}
	stripTimeDiff += (stripTimeMax-stripTimeMin);
	if( currentRecoHit->_trueType==0x01 ) ghostRatio += 1.;
      }
      stripTimeDiff /= (double)tempTrack->_recoHits.size();
      ghostRatio /= (double)tempTrack->_recoHits.size();
      log("debug") << "ghost hit ratio = " << ghostRatio << ", strip time diff = " << stripTimeDiff << endl;
      if( stripTimeDiff>_parGhostTimeDiff ){
	tempTrack->_type += (0x01<<IsGhostTrack);
      }
    }
  }

  //-------------------------//
  // find connectable tracks //
  //-------------------------//
  if( _parConnectTrack ){
    log("debug") << "[Start track connection]" << endl;
    
    Track *targetTrack;
    for(unsigned int itrack=0; itrack<_tempTracks.size(); itrack++){
      tempTrack = _tempTracks.at(itrack);
      if( tempTrack->_type & (0x01<<IsGhostTrack) ){
	log("debug") << "Track " << itrack << " is ghost" << endl;
	continue;
      }
      double lastTime = trackLastTime.at(itrack);
   
      genfit::TrackPoint *tp = fitTrackVec.at(itrack)->getPointWithMeasurementAndFitterInfo(fitTrackVec.at(itrack)->getNumPoints()-1);
      if(tp==NULL){
	continue;
      }

      const int pdg = pdgWanted;
      double chi2optMin=1e6;
      int chi2optMinID = -1;

      for(unsigned int itrack2=0; itrack2<_tempTracks.size(); itrack2++){
	if(itrack==itrack2) continue;
	targetTrack = _tempTracks.at(itrack2);
	if( targetTrack->_type & (0x01<<IsGhostTrack) ) continue;
	double targetTime = targetTrack->_time;
	log("debug") << "current " << itrack << " last time = " << lastTime << ", target " << itrack2 << " start time = " << targetTime << endl;
	if( lastTime<targetTime-_parConnectSearchTime ){
	  log("debug") << "extrapolation is skipped" << endl; 
	  continue;
	}
	
	genfit::TrackPoint *tp2 = fitTrackVec.at(itrack2)->getPointWithMeasurementAndFitterInfo(0);
	if(tp2==NULL){
	  continue;
	}
	genfit::KalmanFittedStateOnPlane kfsop2(*(static_cast<genfit::KalmanFitterInfo*>(tp2->getFitterInfo())->getBackwardUpdate()));

	//-----------------------//
	// forward extrapolation //
	//-----------------------//
	log("debug") << "forward extrapolation" << endl;
	genfit::KalmanFittedStateOnPlane kfsop(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo())->getForwardUpdate()));
	TVector3 posOrg = kfsop.getPos();
	TVector3 momOrg = kfsop.getMom();
	TMatrixDSym covOrg = kfsop.getCov();
	
	genfit::DetPlane *detPlane = 0;
	currentRecoHit = targetTrack->_recoHits.at(0);
	if( currentRecoHit==0 ){
	  log("warning") << "RecoHit is not found in target track" << endl;
	  continue;
	}
	if( !_parUseStripCluster ){
	  detPlane = getRecoHitPlane(currentRecoHit);
	}else{
	  currentStripCluster = currentRecoHit->_stripClusters.at(0);
	  if( currentStripCluster ){
	    detPlane = getStripClusterPlane(currentStripCluster);
	  }else{
	    log("warning") << "StripCluster is not found in current RecoHit" << endl;
	    continue;
	  }
	}
	genfit::SharedPlanePtr plane(detPlane);
	TVector3 uvec = detPlane->getU();
	TVector3 vvec = detPlane->getV();
	TVector3 nvec = detPlane->getNormal();
	TVector3 ovec = detPlane->getO();
	
	// fitted track extrapolation
	genfit::AbsTrackRep *repfit = fitTrackVec.at(itrack)->getCardinalRep();
	repfit->setPropDir(1);
	try{
	  double extrapLenFit = repfit->extrapolateToPlane(kfsop, plane);
	  TMatrixDSym covExp = kfsop.getCov();
	  TVector3 posExp = kfsop.getPos();
	  TVector3 momExp = kfsop.getMom();
	  double uExp = (posExp-ovec)*uvec;
	  double vExp = (posExp-ovec)*vvec;
	  double upExp = momExp*uvec/(momExp*nvec);
	  double vpExp = momExp*vvec/(momExp*nvec);
	  
	  posExp *= 10.; // cm -> mm
	  momExp *= 1000.; // GeV/c -> MeV/c

	  TVector3 posTar = targetTrack->_pos;
	  TVector3 momTar = targetTrack->_p;
	  TMatrixDSym covTar = targetTrack->_covM;
	  log("debug") << "extrap r = " << posExp.Perp() << " +/- " << sqrt(covExp(3,3))*10 << " mm, z = " << posExp.Z() << " +/- " << sqrt(covExp(4,4))*10 << " mm, phi = " << posExp.Phi() << endl;
	  log("debug") << "target r = " << posTar.Perp() << " +/- " << sqrt(covTar(3,3))*10 << " mm, z = " << posTar.Z() << " +/- " << sqrt(covTar(4,4))*10 << " mm, phi = " << posTar.Phi() << endl;
	  log("debug") << "extrap time = " << kfsop.getTime()+tempTrack->_time << " ns" << endl;
	  log("debug") << "target time = " << targetTrack->_time << " ns" << endl;
	  log("debug") << "extrap mom = " << momExp.Mag() << " MeV/c, u' = " << upExp << " +/- " << sqrt(covExp(1,1)) << " v' = " << vpExp << " +/- " << sqrt(covExp(2,2)) << endl;
	  double upTar = momTar*uvec/(momTar*nvec);
	  double vpTar = momTar*vvec/(momTar*nvec);
	  posTar *= 0.1; // mm -> cm
	  double uTar = (posTar-ovec)*uvec;
	  double vTar = (posTar-ovec)*vvec;
	  log("debug") << "target mom = " << momTar.Mag() << " MeV/c, u' = " << upTar << " v' = " << vpTar << endl;
	  
	  double deltaPosRSig = (uTar-uExp)/sqrt(covExp(3,3)+covTar(3,3));
	  double deltaPosZSig = (vTar-vExp)/sqrt(covExp(4,4)+covTar(4,4));
	  double deltaPosSig = deltaPosRSig*deltaPosRSig + deltaPosZSig*deltaPosZSig;
	  double deltaMomSig = (1./momTar.Mag()-1./momExp.Mag())/(sqrt(covExp(0,0)+covTar(0,0))*0.001);
	  double deltaUpSig = (upTar-upExp)/sqrt(covExp(1,1)+covTar(1,1));
	  double deltaVpSig = (vpTar-vpExp)/sqrt(covExp(2,2)+covTar(2,2));
	  double chi2opt = (deltaPosSig + deltaMomSig*deltaMomSig)/3.;
	  double deltaTime = targetTrack->_time - (kfsop.getTime()+tempTrack->_time);
	  log("debug") << "momsig = " << deltaMomSig << " upsig = " << deltaUpSig << " vpsig = " << deltaVpSig << endl;
	  log("debug") << "deltaTime = " << deltaTime << " chi2opt " << chi2opt << endl;
	  
	  if( abs(deltaTime)<_parConnectTimeDiff && chi2opt<chi2optMin ){
	    chi2optMin = chi2opt;
	    chi2optMinID = itrack2;
	  }
	}catch(genfit::Exception& e){
	  log("error") << "Exception in extraplation to plane using genfit::Track" << endl;
	  log("error") << e.what();
	}
	
	//------------------------//
	// backward extrapolation //
	//------------------------//
	log("debug") << "backward extrapolation" << endl;
	
	genfit::DetPlane *detPlane2 = 0;
	currentRecoHit = tempTrack->_recoHits.at(tempTrack->_recoHits.size()-1);
	if( currentRecoHit==0 ){
	  log("warning") << "RecoHit is not found in current track" << endl;
	  continue;
	}
	if( !_parUseStripCluster ){
	  detPlane2 = getRecoHitPlane(currentRecoHit);
	}else{
	  currentStripCluster = currentRecoHit->_stripClusters.at(currentRecoHit->_stripClusters.size()-1);
	  if( currentStripCluster ){
	    detPlane2 = getStripClusterPlane(currentStripCluster);
	  }else{
	    log("warning") << "StripCluster is not found in current RecoHit" << endl;
	    continue;
	  }
	}
	genfit::SharedPlanePtr plane2(detPlane2);
	TVector3 uvec2 = detPlane2->getU();
	TVector3 vvec2 = detPlane2->getV();
	TVector3 nvec2 = detPlane2->getNormal();
	TVector3 ovec2 = detPlane2->getO();

	genfit::AbsTrackRep *repfit2 = fitTrackVec.at(itrack2)->getCardinalRep();
	repfit2->setPropDir(-1);
	try{
	  repfit2->extrapolateToPlane(kfsop2, plane2);
	  TMatrixDSym covExp2 = kfsop2.getCov();
	  TVector3 posExp2 = kfsop2.getPos();
	  TVector3 momExp2 = kfsop2.getMom();
	  double uExp2 = (posExp2-ovec2)*uvec2;
	  double vExp2 = (posExp2-ovec2)*vvec2;
	  double upExp2 = momExp2*uvec2/(momExp2*nvec2);
	  double vpExp2 = momExp2*vvec2/(momExp2*nvec2);

	  double uOrg = (posOrg-ovec2)*uvec2;
	  double vOrg = (posOrg-ovec2)*vvec2;
	  double upOrg = momOrg*uvec2/(momOrg*nvec2);
	  double vpOrg = momOrg*vvec2/(momOrg*nvec2);
	  posExp2 *= 10.; // cm -> mm
	  posOrg *= 10.;
	  momExp2 *= 1000.; // GeV/c -> MeV/c
	  momOrg *= 1000.;
	  log("debug") << "extrap r = " << posExp2.Perp() << " +/- " << sqrt(covExp2(3,3))*10 << " mm, z = " << posExp2.Z() << " +/- " << sqrt(covExp2(4,4))*10 << " mm, phi = " << posExp2.Phi() << endl;
	  log("debug") << "target r = " << posOrg.Perp() << " +/- " << sqrt(covOrg(3,3))*10 << " mm, z = " << posOrg.Z() << " +/- " << sqrt(covOrg(4,4))*10 << " mm, phi = " << posOrg.Phi() << endl;
	  log("debug") << "extrap time = " << targetTrack->_time+kfsop2.getTime() << " ns" << endl;
	  double deltaPosRSig2 = (uOrg-uExp2)/sqrt(covOrg(3,3)+covExp2(3,3));
	  double deltaPosZSig2 = (vOrg-vExp2)/sqrt(covOrg(4,4)+covExp2(4,4));
	  double deltaPosSig2 = deltaPosRSig2*deltaPosRSig2 + deltaPosZSig2*deltaPosZSig2;
	  double deltaMomSig2 = (1./momOrg.Mag()-1./momExp2.Mag())/(sqrt(covOrg(0,0)+covExp2(0,0))*0.001);
	  double deltaUpSig2 = (upOrg-upExp2)/sqrt(covOrg(1,1)+covExp2(1,1));
	  double deltaVpSig2 = (vpOrg-vpExp2)/sqrt(covOrg(2,2)+covExp2(2,2));
	  double deltaTimeBack = lastTime-(targetTrack->_time+kfsop2.getTime());
	  double chi2optBack = (deltaPosSig2 + deltaMomSig2*deltaMomSig2)/3.;
	  log("debug") << "deltaTimeBack = " << deltaTimeBack << ", chi2optBack " << chi2optBack << endl;
	  if( abs(deltaTimeBack)<_parConnectTimeDiff && chi2optBack<chi2optMin ){
	    chi2optMin = chi2optBack;
	    chi2optMinID = itrack2;
	  }
	}catch(genfit::Exception& e){
	  log("error") << "Exception in backward extraplation to plane using genfit::Track" << endl;
	  log("error") << e.what();
	}
      }
    
      if( chi2optMin<_parConnectChi2NDFMax ){
	log("debug") << "End of track " << itrack << " is connected to track " << chi2optMinID << endl;
	tempTrack->_next = _tempTracks.at(chi2optMinID);
	_tempTracks.at(chi2optMinID)->_prev = tempTrack;
      }
    }

  }

  for(unsigned int itrack=0; itrack<_tempTracks.size(); itrack++){
    _fittedTracks.push_back(_tempTracks.at(itrack));
  }

  for(unsigned int itrack=0; itrack<fitTrackVec.size(); itrack++){
    delete fitTrackVec.at(itrack);
  }
  fitTrackVec.clear();
}


void TrkFitting::finish(Parameter *)
{
  if(_fitter) delete _fitter;

  if(_parShowEventDisplay && _display){
    _display->setOptions("ABDEFGHMPT");
    _display->open();
  }
}

void TrkFitting::setRecoHitPosition(genfit::Track *fitTrack, const RecoHit *recoHit)
{
  const double resolution = _parStripPitch*0.1/sqrt(12); // cm
  TMatrixDSym hitCov(2);
  hitCov.UnitMatrix();
  hitCov *= (resolution*resolution);
  
  double phi = recoHit->_pos.Phi();
  TVectorD hitCoords(2);
  hitCoords[0] = recoHit->_pos.Perp()*0.1; // mm -> cm
  hitCoords[1] = recoHit->_pos.Z()*0.1; // mm -> cm
  int sensorID = 0;
  if( recoHit->_stripClusters.size()>0 && recoHit->_stripClusters.at(0) ){
    sensorID = recoHit->_stripClusters.at(0)->sensorID();
  }
    
  int hitId = fitTrack->getNumPoints();
  genfit::PlanarMeasurement *measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, sensorID, hitId, nullptr );
  TVector3 o(0.,0.,0.);
  TVector3 u(cos(phi), sin(phi), 0.);
  TVector3 v(0., 0., 1.);
  measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u,v, new genfit::MyFinitePlane())), 0 );
  fitTrack->insertPoint(new genfit::TrackPoint(measurement,fitTrack));
}

void TrkFitting::setStripClusterPosition(genfit::Track *fitTrack, const StripCluster *stripCluster)
{
  double originX, originY, originZ;
  unsigned int vaneID = stripCluster->vaneID();
  double phi = -2*TMath::Pi()*vaneID/(double)_parNVane;
  unsigned int sensorID = stripCluster->sensorID();

  const double resolution = _parStripPitch*0.1/sqrt(12); // cm
  TMatrixDSym hitCov(1);
  hitCov.UnitMatrix();
  hitCov *= (resolution*resolution);
  int hitId = fitTrack->getNumPoints();

  if( sensorID&0x01 ){
    originX = _parSensorOriginX.at(0) * 0.1; // cm
  }else{
    originX = _parSensorOriginX.at(1) * 0.1; // cm
  }
  if( sensorID&0x02 ){
    originZ = _parSensorOriginZ.at(0) * 0.1; // cm
  }else{
    originZ = _parSensorOriginZ.at(1) * 0.1; // cm
  }
  if( sensorID&0x04 ){
    originY = _parSensorDistance*0.5 * 0.1; // cm
  }else{
    originY = -_parSensorDistance*0.5 * 0.1; // cm
  }

  unsigned int stripID = stripCluster->_stripID;
  TVectorD hitCoords(1);
  hitCoords[0] = ((_parNStripPerRow-1-stripCluster->_pos)* _parStripPitch + _parPosFirstStrip) * 0.1; // cm
	    
  if( sensorID&0x08 ){
    originZ = -originZ;
    if( stripCluster->isZStrip() ){
      hitCoords[0] = -hitCoords[0];
    }
  }	    
  genfit::PlanarMeasurement *measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, (int)sensorID, hitId, nullptr );       
  if( stripCluster->isZStrip() ){
    measurement->setStripV();
  }
  
  TVector3 o(originX*cos(phi)-originY*sin(phi), originX*sin(phi)+originY*cos(phi), originZ);  
  TVector3 u(cos(phi), sin(phi), 0.);
  TVector3 v(0., 0., 1.);
  measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u,v, new genfit::MyFinitePlane()) ), 0 );      
  fitTrack->insertPoint(new genfit::TrackPoint(measurement,fitTrack) );
}


genfit::DetPlane* TrkFitting::getRecoHitPlane(const RecoHit *recoHit)
{
  double phi = recoHit->_pos.Phi();
  TVector3 o(0., 0., 0.);
  TVector3 u(cos(phi), sin(phi), 0.);
  TVector3 v(0., 0., 1.);
  return new genfit::DetPlane(o, u, v, new genfit::MyFinitePlane());
}


genfit::DetPlane* TrkFitting::getStripClusterPlane(const StripCluster *stripCluster)
{
  double originX, originY, originZ;
  unsigned int vaneID = stripCluster->vaneID();
  double phi = -2*TMath::Pi()*vaneID/(double)_parNVane;
  unsigned int sensorID = stripCluster->sensorID();

  if( sensorID&0x01 ){
    originX = _parSensorOriginX.at(0) * 0.1; // cm
  }else{
    originX = _parSensorOriginX.at(1) * 0.1; // cm
  }
  if( sensorID&0x02 ){
    originZ = _parSensorOriginZ.at(0) * 0.1; // cm
  }else{
    originZ = _parSensorOriginZ.at(1) * 0.1; // cm
  }
  if( sensorID&0x04 ){
    originY = _parSensorDistance*0.5 * 0.1; // cm
  }else{
    originY = -_parSensorDistance*0.5 * 0.1; // cm
  }

  if( sensorID&0x08 ){
    originZ = -originZ;
  }
  
  TVector3 o(originX*cos(phi)-originY*sin(phi), originX*sin(phi)+originY*cos(phi), originZ);  
  TVector3 u(cos(phi), sin(phi), 0.);
  TVector3 v(0., 0., 1.);
  return new genfit::DetPlane(o,u,v, new genfit::MyFinitePlane());
}

double TrkFitting::getTrueTime(const RecoHit *recoHit, const bool doFirst)
{
  double trueTimeMin = 1e9;
  double trueTimeMax = -1e9;
  for(unsigned int iClus=0; iClus<recoHit->_stripClusters.size(); iClus++){
    double trueTime = getTrueTime(recoHit->_stripClusters.at(iClus), doFirst );
    if( trueTime<trueTimeMin ){
      trueTimeMin = trueTime;
    }
    if( trueTime>trueTimeMax ){
      trueTimeMax = trueTime;
    }
  }
  if( doFirst ) return trueTimeMin;
  else return trueTimeMax;
}

double TrkFitting::getTrueTime(const StripCluster *stripCluster, const bool doFirst)
{
  double trueTimeMin = 1e9;
  double trueTimeMax = -1e9;

  for(unsigned int iStrip=0; iStrip<stripCluster->_stripHits.size(); iStrip++){	      
    for(unsigned int iSim=0; iSim<stripCluster->_stripHits.at(iStrip)->_simHits.size(); iSim++){
      const SimHit *simhit = stripCluster->_stripHits.at(iStrip)->_simHits.at(iSim);
      if(simhit){
	double trueTime = simhit->_time;
	if( trueTime<trueTimeMin ){
	  trueTimeMin = trueTime;
	}
	if( trueTime>trueTimeMax ){
	  trueTimeMax = trueTime;
	}
      }
    }
  }
  if( doFirst ) return trueTimeMin;
  else return trueTimeMax;
}
