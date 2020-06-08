// RunManager.cxx

// TreeProc includes
#include "TreeProc/RunManager.h"
#include "TreeProc/ProcManager.h"
#include "TreeProc/TreeManager.h"
#include "TreeProc/LogStream.h"

#include <iostream>

#include "TSystem.h"

using namespace std;
using namespace TreeProc;

RunManager * RunManager::_instance = 0;

RunManager * RunManager::instance(){
  if(_instance) return _instance;
  else return _instance = new RunManager;
}

RunManager::RunManager() : _arranger (0), _evIn(0), _evOut(0)
{}

void RunManager::loadLibraries(std::vector<std::string> &libs)
{
  for(unsigned int i=0;i<libs.size();i++){
    gSystem->Load(libs[i].c_str());
  }
}

void RunManager::config(Parameter *param, const char *arranger){
  _param = param;
  _parArrangerName = arranger;
  
  // _arranger
  if(_parArrangerName){
    _arranger = ObjectFactory<EventArranger>::instance()->makeObject(_parArrangerName, _parArrangerName);
    if(_arranger){
      log() << "EventArranger " << _parArrangerName << " has loaded." << endl;
      Parameter *paramArr = ProcManager::instance()->getParam(_parArrangerName, _parArrangerName);
      _arranger->config(paramArr);
    }
    else
      log() << "EventArranger " << _parArrangerName << " has been failed to load. Using default arrangement." << endl;
  }
}

int RunManager::getNextEntry(){
  bool b = TreeProc::TreeManager::instance()->getEntry(_evIn++);
  if(!b)_evIn = -1;
  return _evIn;
}

void RunManager::eventLoop(){
  // loading global configurations
  _parEventNumToProcess = _param->get("EventNumToProcess", int(-1));
  _parEventIdxToStart = _param->get("EventIdxToStart", int(-1));
  _parWritePeriod = _param->get("WritePeriod", int(100));
  _parPrintPeriod = _param->get("PrintPeriod", int(100));
  _parInputFile = _param->get("InputFile", string(""));
  
  int nStop = 0;
  if (_parInputFile != ""){
    int nEv = TreeManager::instance()->getEntries();
    int nEvent = _parEventNumToProcess;
    if(nEv < nEvent || nEvent == -1) nEvent = nEv;
    _evIn = _parEventIdxToStart;
    nStop = _evIn + nEvent;
  }
  else {
    _evIn = 0;
    nStop = _parEventNumToProcess;
  }
  
  log() << "RunManager::eventLoop(): event " << _evIn << " to " << nStop - 1 << " will be processed" << endl;
  
  int prev_evIn = 0;
  while(_evIn < nStop){
    if(prev_evIn == 0 || _evIn - prev_evIn > _parPrintPeriod){
      log() << "Event " << _evIn << " processing..." << endl;
      prev_evIn = _evIn;
    }

    if(_arranger){
      Parameter *param = ProcManager::instance()->getParam(_parArrangerName, _parArrangerName);
      bool b = _arranger->arrange(_evIn, param);
      if(!b)continue;
    }else{
      getNextEntry();
      if(_evIn == -1)break;
    }
    ProcManager::instance()->process(_evOut);
    TreeManager::instance()->fill();

    if(_evOut % _parWritePeriod == _parWritePeriod - 1)
      TreeManager::instance()->write();
    
    _evOut ++;
  }
  
  TreeManager::instance()->write();

  // calling config: todo: consider if this is suitable here
  ProcManager::instance()->finish();
}
