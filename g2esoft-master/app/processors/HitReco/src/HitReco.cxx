// HitReco.cxx

#include "TreeProc/TreeManager.h"
#include "TreeProc/LogStream.h"
#include "HitReco.h"
//#include "g2track/g2track.h"

#include <iostream>

using namespace g2esoft;
using namespace TreeProc;
using namespace std;

static Factory<Processor, HitReco> _f("HitReco");

HitReco::HitReco(const char *procName, const char *instName) : Processor(procName, instName)
{
}

HitReco::~HitReco(){}

void HitReco::config(Parameter *param){
  _parRecoHitName = param->get("RecoHitName",string("RecoHits"));
  TreeManager::instance()->registerBranch(_parRecoHitName.c_str(), &_recoHits);
  _parStripClusterName = param->get("StripClusterName",string("StripClusters"));

  _parNVane = param->get("NVane",(int)40);
  _parDoubleRow = param->get("DoubleRow",(bool)true);
  _parSensorOriginX = param->get("SensorOriginX",vector<double>{90.96,190.23});
  _parSensorOriginZ = param->get("SensorOriginZ",vector<double>{0.5,99.77});
  _parStripPitch = param->get("StripPitch",(double)0.190);
  _parSensorSize = param->get("SensorSize",(double)98.77);
  _parInactiveMiddle = param->get("InactiveMiddle",(double)0.5);
  _parTimeDiv = param->get("TimeDiv",(int)5);
  _parTimeMatch = param->get("TimeMatch", (double)5.);

}

void HitReco::process(int nEvent, Parameter *){
  log("debug") << "start HitReco process" << endl;

  // input: stripClusters
  const vector<const StripCluster *> &stripClusters = TreeManager::instance()->getBranchVec<StripCluster>(_parStripClusterName.c_str());

  // clear recoHits
  for(unsigned int i=0; i<_recoHits.size(); i++){
    delete _recoHits[i];
  }
  _recoHits.clear();

  std::vector<int> stripClusterZIndices;
  std::vector<int> stripClusterRIndices;

  for(unsigned int i=0; i<stripClusters.size(); i++){
    if( this->isZStrip( stripClusters.at(i)->_stripID ) ){
      stripClusterZIndices.push_back(i);
    }else{
      stripClusterRIndices.push_back(i);
    }
  }

  TVector3 pos;

  for(unsigned int i=0; i<stripClusterZIndices.size(); i++){
    int indexZ = stripClusterZIndices.at(i);
    // strip ID to 3D position
    unsigned int stripIDZ = stripClusters.at(indexZ)->_stripID;
    unsigned int vaneIDZ = this->vaneID( stripIDZ );

    // sensor ID without Z or R
    unsigned int sensorIDZwo = (this->sensorID( stripIDZ ))&0x0b;
    
    // strip only ID
    unsigned int stripOnlyIDZ = this->stripOnlyID( stripIDZ );
    
    // strip only ID for GSB 2 bit
    unsigned int stripOnlyGSBIDZ = stripOnlyIDZ/0x100;

    pos = this->stripPosition( stripIDZ, stripClusters.at(indexZ)->_pos );
    double posZ = pos.Z();
    double posPhi = pos.Phi();

    log("debug") << "Z strip: vaneID = " << vaneIDZ << " sensorIDZwo = " << sensorIDZwo << " stripOnlyGSBID = " << stripOnlyGSBIDZ << " time = " << stripClusters.at(indexZ)->_time << endl;

    for(unsigned int j=0; j<stripClusterRIndices.size(); j++){
      int indexR = stripClusterRIndices.at(j);

      unsigned int stripIDR = stripClusters.at(indexR)->_stripID;
      unsigned int vaneIDR = this->vaneID( stripIDR );
      unsigned int sensorIDRwo = ((this->sensorID( stripIDR )) & 0x0b);
      unsigned int stripOnlyIDR = this->stripOnlyID( stripIDR );
      unsigned int stripOnlyGSBIDR = stripOnlyIDR/0x100;
      double posR = this->stripPosition( stripIDR, stripClusters.at(indexR)->_pos ).Perp();

      log("debug") << "R strip: vaneID = " << vaneIDR << " sensorIDRwo = " << sensorIDRwo << " stripOnlyGSBID = " << stripOnlyGSBIDR << " time = " << stripClusters.at(indexR)->_time << endl;

      if( abs((int)(stripClusters.at(indexZ)->_time) - (int)(stripClusters.at(indexR)->_time))<=_parTimeMatch && // hits are regarded as same origin
	  vaneIDZ==vaneIDR && sensorIDZwo==sensorIDRwo &&
	  (!_parDoubleRow ||
	   (((stripOnlyGSBIDZ==0 && stripOnlyGSBIDR==0) ||
	     (stripOnlyGSBIDZ==1 && stripOnlyGSBIDR==2) ||
	     (stripOnlyGSBIDZ==2 && stripOnlyGSBIDR==1) ||
	     (stripOnlyGSBIDZ==3 && stripOnlyGSBIDR==3)))) ){

	log("debug") << "matched hits" << endl;

	RecoHit *recoHit = new RecoHit;
	recoHit->_stripClusters.push_back( stripClusters.at(indexZ) );
	recoHit->_stripClusters.push_back( stripClusters.at(indexR) );
	recoHit->_pos.SetXYZ(posR*cos(posPhi), posR*sin(posPhi), posZ);
	recoHit->_time = 0.5*(stripClusters.at(indexZ)->_time + stripClusters.at(indexR)->_time);
	//recoHit->_tot = 0.5*(stripClusters.at(indexZ)->_tot + stripClusters.at(indexR)->_tot);
	// how to define?
	//recoHit->isNoize = false;

	log("debug") << "recoHit pos(r,phi,z) = " << recoHit->_pos.Perp() << ", " << recoHit->_pos.Phi() << ", " << recoHit->_pos.Z() << endl;

	bool isGhost = true;
	const MCParticle *mcp1 = 0;
	const MCParticle *mcp2 = 0;
	double time1, time2;
	for(unsigned int k1=0; k1<stripClusters.at(indexZ)->_stripHits.size(); ++k1){
	  for(unsigned int k2=0; k2<stripClusters.at(indexZ)->_stripHits.at(k1)->_simHits.size(); ++k2){
	    if( stripClusters.at(indexZ)->_stripHits.at(k1)->_simHits.at(k2)->getMCStep() ){
	      mcp1 = stripClusters.at(indexZ)->_stripHits.at(k1)->_simHits.at(k2)->getMCStep()->getMCP();
	      time1 = stripClusters.at(indexZ)->_stripHits.at(k1)->_simHits.at(k2)->_time;
	    }else{
	      continue;
	    }

	    for(unsigned int l1=0; l1<stripClusters.at(indexR)->_stripHits.size(); ++l1){
	      for(unsigned int l2=0; l2<stripClusters.at(indexR)->_stripHits.at(l1)->_simHits.size(); ++l2){
		if( stripClusters.at(indexR)->_stripHits.at(l1)->_simHits.at(l2)->getMCStep() ){
		  mcp2 = stripClusters.at(indexR)->_stripHits.at(l1)->_simHits.at(l2)->getMCStep()->getMCP();
		  time2 = stripClusters.at(indexR)->_stripHits.at(l1)->_simHits.at(l2)->_time;
		}else{
		  continue;
		}
		
		if( mcp1 && mcp1==mcp2 && abs(time1-time2)<0.1 ){
		  isGhost = false;
		  break;
		}
	      }
	      if( !isGhost ) break;
	    }
	    if( !isGhost ) break;
	  }
	  if( !isGhost ) break;
	}

	if( isGhost ){
	  recoHit->_trueType |= (0x01<<IsGhostHit);
	}

	log("debug") << "isGhost = " << (recoHit->_trueType&(0x01<<IsGhostHit)) << endl;

	_recoHits.push_back( recoHit );
      }
    }
  }

  log("info") << _recoHits.size() << " reco hits created." << endl; 
}

void HitReco::finish(Parameter *){}


inline bool HitReco::isZStrip(const unsigned int stripID)
{
  return ((this->sensorID(stripID)) & 0x04);
}

inline bool HitReco::isRStrip(const unsigned int stripID)
{
  return !(this->isZStrip(stripID));
}

inline unsigned int HitReco::vaneID(const unsigned int stripID)
{
  return (stripID/0x1000000);
}

inline unsigned int HitReco::sensorID(const unsigned int stripID)
{
  return (stripID%0x1000000)/0x10000;
}

inline unsigned int HitReco::stripOnlyID(const unsigned int stripID)
{
  return (stripID%0x10000);
}

TVector3 HitReco::stripPosition(const unsigned int stripID)
{
  unsigned int vaneID = this->vaneID(stripID);

  double phi = -2*TMath::Pi()*vaneID/(double)_parNVane; // ranges from -2pi to 0

  unsigned int sensorID = this->sensorID(stripID);
  double sensorR = _parSensorOriginX.at((sensorID&0x01)?0:1) + 0.5*_parSensorSize;
  double sensorZ = _parSensorOriginZ.at((sensorID&0x02)?0:1) + 0.5*_parSensorSize;
  
  unsigned int stripOnlyID = this->stripOnlyID( stripID );
  double length = 0.5*(_parInactiveMiddle+0.5*_parSensorSize);
  if( this->isRStrip(stripID) ){
    sensorR += ((255.5-(double)(stripOnlyID&0x1ff))*_parStripPitch);
    if( (stripOnlyID&0x200)==0 ){
      sensorZ += length;
    }else{
      sensorZ -= length;
    }
  }else{
    if( (stripOnlyID&0x200)==0 ){
      sensorR += length;
    }else{
      sensorR -= length;
    }
    sensorZ += ((255.5-(double)(stripOnlyID&0x1ff))*_parStripPitch);
  }

  if( this->sensorID( stripID )&0x08 ) sensorZ = -sensorZ;

  TVector3 pos(sensorR*cos(phi), sensorR*sin(phi), sensorZ);
  return pos;
}


TVector3 HitReco::stripPosition(const unsigned int stripID, const double stripPos)
{
  unsigned int vaneID = this->vaneID(stripID);

  double phi = -2*TMath::Pi()*vaneID/(double)_parNVane; // ranges from -2pi to 0

  unsigned int sensorID = this->sensorID(stripID);
  double sensorR = _parSensorOriginX.at((sensorID&0x01)?0:1) + 0.5*_parSensorSize;
  double sensorZ = _parSensorOriginZ.at((sensorID&0x02)?0:1) + 0.5*_parSensorSize;

  unsigned int stripOnlyID = this->stripOnlyID( stripID );
  double length = 0.5*(_parInactiveMiddle+0.5*_parSensorSize);
  if( this->isRStrip(stripID) ){
    sensorR += ((255.5-stripPos)*_parStripPitch);
    if( (stripOnlyID&0x200)==0 ){
      sensorZ += length;
    }else{
      sensorZ -= length;
    }
  }else{
    if( (stripOnlyID&0x200)==0 ){
      sensorR += length;
    }else{
      sensorR -= length;
    }
    sensorZ += ((255.5-stripPos)*_parStripPitch);
  }

  if( this->sensorID( stripID )&0x08 ) sensorZ = -sensorZ;

  TVector3 pos(sensorR*cos(phi), sensorR*sin(phi), sensorZ);
  return pos;
}
