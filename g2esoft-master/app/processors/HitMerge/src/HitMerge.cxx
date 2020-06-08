// HitMerge.cxx

#include "TreeProc/TreeManager.h"
#include "TreeProc/LogStream.h"
#include "HitMerge.h"
//#include "g2track/g2track.h"

#include <iostream>
#include <map>
#include <unordered_map>
#include <set>

using namespace g2esoft;
using namespace TreeProc;
using namespace std;

static Factory<Processor, HitMerge> _f("HitMerge");

HitMerge::HitMerge(const char *procName, const char *instName) : Processor(procName, instName)
{}

HitMerge::~HitMerge(){}

void HitMerge::config(Parameter *param)
{
  _parStripClusterName = param->get("StripClusterName",string("StripClusters"));
  TreeManager::instance()->registerBranch(_parStripClusterName.c_str(),&_stripClusters);
  _parStripHitName = param->get("StripHitName",string("StripHits"));
  _parTimeMatch = param->get("TimeMatch",(double)1.);
}

void HitMerge::process(int nEvent, Parameter *param)
{
  log("debug") << "Come to event " << nEvent << endl;

  const vector<const StripHit *> &stripHits = TreeManager::instance()->getBranchVec<StripHit>(_parStripHitName.c_str());

  for(unsigned int i=0; i<_stripClusters.size(); i++){
    delete _stripClusters[i];
  }
  _stripClusters.clear();

  multimap<unsigned int, const StripHit*> mapStripHits;
  typedef multimap<unsigned int, const StripHit*>::value_type mapvaltyp;
  set<unsigned int> keyset;

  for(unsigned int i=0; i<stripHits.size(); ++i){
    unsigned int stripID = stripHits.at(i)->_stripID;
    stripID &= 0xffffffe00;
    mapStripHits.insert(mapvaltyp(stripID, stripHits.at(i)));
    keyset.insert(stripID);
  }

  typedef unordered_map<const StripHit *, StripCluster *>::value_type uomaptyp;

  for(auto keyit=keyset.begin(); keyit!=keyset.end(); ++keyit){
    auto p = mapStripHits.equal_range(*keyit);
    unordered_map<const StripHit *, StripCluster *> mapStripHitCluster;
    vector<StripCluster *> vecStripCluster;
    for(auto it=p.first; it!=p.second; ++it){
      if( mapStripHitCluster.find(it->second)==mapStripHitCluster.end() ){
	StripCluster *stripCluster = new StripCluster;
	stripCluster->_stripID = *keyit;
	stripCluster->_stripHits.push_back( it->second );
	mapStripHitCluster.insert(uomaptyp(it->second, stripCluster));
	vecStripCluster.push_back( stripCluster );
      }
      auto it2 = it;
      it2++;
      for(; it2!=p.second; ++it2){
	if( mapStripHitCluster.find(it2->second)==mapStripHitCluster.end() && 
	    abs(it->second->_time - it2->second->_time)<=_parTimeMatch &&
	    abs((int)it->second->_stripID-(int)it2->second->_stripID)<=1 ){
	  mapStripHitCluster.insert(uomaptyp(it2->second, mapStripHitCluster.at(it->second)));
	  mapStripHitCluster.at(it->second)->_stripHits.push_back( it2->second );
	}
      }
    }
  
    for(auto it=vecStripCluster.begin(); it!=vecStripCluster.end(); ++it){
      double posAve = 0.;
      double timeAve = 0.;
      double edepSum = 0.;
      for(unsigned int i=0; i<(*it)->_stripHits.size(); ++i){
	double edep = (*it)->_stripHits.at(i)->_tot/100.*0.087; // MeV 
	posAve += ((*it)->_stripHits.at(i)->_stripID & 0x1ff)*edep;
	timeAve += (*it)->_stripHits.at(i)->_time * edep;
	edepSum += edep;
      }
      if( edepSum>0. ){
	posAve /= edepSum;
	timeAve /= edepSum;
      }
      (*it)->_pos = posAve;
      (*it)->_time = timeAve;
      (*it)->_edep = edepSum;
      if( posAve>255.5 ){
	(*it)->_stripID += 0x100;
      }

      _stripClusters.push_back( *it );
    }
  }
  
  log("info") << _stripClusters.size() << " strip clusters created" << endl;
}

void HitMerge::finish(Parameter *){}

