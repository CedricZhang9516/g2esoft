#include "G4SimSteppingAction.h"

#include "TreeProc/LogStream.h"

#include "G4VProcess.hh"
#include "G4SteppingManager.hh"

using namespace TreeProc;
using namespace g2esoft;
using namespace std;

G4SimSteppingAction::G4SimSteppingAction (std::vector<MCParticle *> &mcParticles, std::vector<const MCStep *> &mcSteps, std::vector<const SimHit *> &simHits) 
  : _mcParticles (mcParticles), 
    _mcSteps (mcSteps), 
    _simHits (simHits), 
    _currentEventID(0), 
    _currentTrackID(0), 
    _currentPdg(0), 
    _edepSum(0.), 
    _preTime(0.) {
  log("debug") << "\x1b[33m[G4Sim] G4SimSteppingAction constructor starts.\x1b[0m" << endl;
  //log("debug") << "\x1b[33m[G4Sim] G4SimSteppingAction constructor done.\x1b[0m" << endl;
}

G4SimSteppingAction::~G4SimSteppingAction () {
  log("debug") << "\x1b[33m[G4Sim] G4SimSteppingAction destructor starts.\x1b[0m" << endl;
  //log("debug") << "\x1b[33m[G4Sim] G4SimSteppingAction destructor done.\x1b[0m" << endl;
}

void G4SimSteppingAction::UserSteppingAction (const G4Step *step) {
  const G4String processName = step->GetPostStepPoint ()->GetProcessDefinedStep ()->GetProcessName ();
  G4Track *track = step->GetTrack ();
  const G4double gtime = track->GetGlobalTime(); // ns

  if (gtime>40000.) {
    log("info") << "GlobalTime: " << gtime << ", exceeds 40 us" << endl;
    track->SetTrackStatus(fStopAndKill);
  }

  if (processName == "DecayWithSpin") {
    const G4TrackVector *secondary = fpSteppingManager->GetSecondary ();
    const unsigned int lastID = _mcParticles.size();
    const G4int trackID = track->GetTrackID();
    
    if(lastID>0){
      _mcpMap[trackID] = _mcParticles.at(lastID-1);
    }
    
    G4int trackID_sec = trackID;
    for (G4TrackVector::const_iterator p = secondary->begin (); p != secondary->end (); ++p) {
      const G4ThreeVector x_sec = (*p)->GetPosition ();
      const G4ThreeVector p_sec = (*p)->GetMomentum ();
      const G4ThreeVector s_sec = (*p)->GetPolarization ();
      trackID_sec++;
      MCParticle *mcp_secondary = new MCParticle ();
      mcp_secondary->_pdg = (*p)->GetDefinition ()->GetPDGEncoding ();
      mcp_secondary->_p = TLorentzVector (p_sec.x (), p_sec.y (), p_sec.z (), (*p)->GetTotalEnergy ());
      mcp_secondary->_pol = TVector3 (s_sec.x (), s_sec.y (), s_sec.z ());
      mcp_secondary->_time = (*p)->GetGlobalTime ();
      mcp_secondary->_prodVertex = TVector3 (x_sec.x (), x_sec.y (), x_sec.z ());
      mcp_secondary->_muonID = _currentEventID;
      if(lastID>0){
	mcp_secondary->_parent = _mcParticles.at(lastID-1);
	_mcParticles.at(lastID-1)->_daughters.push_back( mcp_secondary );
      }
      _mcpMap[trackID_sec] = mcp_secondary;
      _mcParticles.push_back (mcp_secondary);
    }
  }
  else {
    const G4StepPoint *preStep = step->GetPreStepPoint();
    const G4VPhysicalVolume *physVolume = preStep->GetPhysicalVolume();
    const G4String volumeName = physVolume->GetName();

    if(volumeName=="Sensor"){
      const G4StepPoint *postStep = step->GetPostStepPoint();      
      const G4StepStatus pointIn = preStep->GetStepStatus();
      const G4StepStatus pointOut = postStep->GetStepStatus();
      const G4double edep = step->GetTotalEnergyDeposit();
      const G4int trackID = track->GetTrackID ();

      switch(pointIn){
      case fGeomBoundary:
      case fAlongStepDoItProc:
      case fUndefined:
	if( _currentTrackID!=0 ){
	  log("debug") << "trackID was not reset. add previous hit" << std::endl;
	  if( _edepSum>0. ){
	    CreateHit();
	  }
	}
	_edepSum = edep;
	_currentTrackID = trackID;
	_currentPdg = track->GetDefinition ()->GetPDGEncoding ();
	_prePosition = preStep->GetPosition();
	_preTime = preStep->GetGlobalTime();
	_preMomentum = preStep->GetMomentum();
	_preEnergy = preStep->GetTotalEnergy();
	_prePolarization = preStep->GetPolarization();
	break;
      case fPostStepDoItProc:
      default:
	_edepSum += edep;
	break;
      }

      _postPosition = postStep->GetPosition();
      
      switch(pointOut){
      case fGeomBoundary: 
      case fAlongStepDoItProc:{
	if( _currentTrackID!=trackID ) log("warnng") << "trackID is  different" << std::endl;
	if( _edepSum>0. ){
	  CreateHit();
	}
	_currentTrackID = 0;
	_edepSum = 0.;
	break;
      }
      case fPostStepDoItProc:
      default:
	break;
      }

    }
  }
}

void G4SimSteppingAction::ClearEvent()
{
  _currentEventID = 0;
  _currentTrackID = 0;
  _currentPdg = 0;
  _edepSum = 0.;
  for(auto it_map=_mcpMap.begin(); it_map!=_mcpMap.end(); it_map++){
    delete it_map->second;
  }
  _mcpMap.clear();
}

void G4SimSteppingAction::CreateHit()
{
  SimHit *hit = new SimHit ();
  hit->_pos[0] = TVector3 (_prePosition.x(), _prePosition.y(), _prePosition.z());
  hit->_pos[1] = TVector3 (_postPosition.x (), _postPosition.y (), _postPosition.z ());
  hit->_time = _preTime;
  hit->_edep = _edepSum;

  MCStep *mcStep = new MCStep ();
  mcStep->_pos = hit->_pos[0];
  mcStep->_time = hit->_time;
  mcStep->_p = TLorentzVector(_preMomentum.x(), _preMomentum.y(), _preMomentum.z(), _preEnergy);
  auto it_map = _mcpMap.find(_currentTrackID);
  if(it_map != _mcpMap.end()){
    mcStep->_mcp = it_map->second;
    it_map->second->_steps.push_back( mcStep );
  }else{
    MCParticle *mcp = new MCParticle ();
    mcp->_pdg = _currentPdg;
    mcp->_p = mcStep->_p;
    mcp->_pol = TVector3 (_prePolarization.x(), _prePolarization.y(), _prePolarization.z());
    mcp->_time = mcStep->_time;
    mcp->_prodVertex = mcStep->_pos;
    mcp->_muonID = _currentEventID;
    mcp->_steps.push_back( mcStep );
    _mcParticles.push_back(mcp);
    _mcpMap[_currentTrackID] = mcp;
    mcStep->_mcp = mcp;
  }
  hit->_mcStep = mcStep;
  _simHits.push_back (hit);
  _mcSteps.push_back (mcStep);
}
