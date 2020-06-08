// ProcManager.cxx

// ROOT includes
#include "TApplication.h"

// TreeProc includes
#include "TreeProc/ProcManager.h"
#include "TreeProc/Processor.h"
#include "TreeProc/Parameter.h"
#include "TreeProc/LogStream.h"

#include <iostream>
#include <stdexcept>

using namespace std;
using namespace TreeProc;

// static var/func for singleton
ProcManager * ProcManager::_instance = 0;
ProcManager * ProcManager::instance(){
  if(_instance) return _instance;
  else return _instance = new ProcManager;
}

ProcManager::ProcManager()
{}

void ProcManager::initializeProc(const char *procName, const char *instanceName){
  if(_procMap.find(instanceName) != _procMap.end()){
    throw runtime_error("ProcManager::initializeProc(): process with the same name already exists.");
  }

  Processor *proc = ObjectFactory<Processor>::instance()->makeObject(procName, instanceName);
  _procMap[instanceName] = proc;
  _processors.push_back(procPair(instanceName, proc));
}

// note: for_each is not easily applicable for map (returns pairs)

void ProcManager::config(Parameter *rootParams)
{
  _rootParams = rootParams;
  
  procVecType::iterator it;
  for(it = _processors.begin(); it != _processors.end(); it++){
    Processor *proc = it->second;
    Parameter *param = getParam(proc->getProcName(), proc->getInstName());
    log().config().setProc(proc->getProcName(),proc->getInstName());
    proc->config(param);
    log().config().disableProc();
  }
}

void ProcManager::process(int nEvent)
{
  //cout << "Process event # " << nEvent << endl;
  //cout << "map size = " << _procMap.size() << endl;
  procVecType::iterator it;
  for(it = _processors.begin(); it != _processors.end(); it++){
    Processor *proc = it->second;
    Parameter *param = getParam(proc->getProcName(), proc->getInstName());

    log().config().setProc(proc->getProcName(),proc->getInstName());
    proc->process(nEvent, param);
    log().config().disableProc();
  }
}

void ProcManager::finish()
{
  procVecType::iterator it;
  for(it = _processors.begin(); it != _processors.end(); it++){
    Processor *proc = it->second;
    Parameter *param = getParam(proc->getProcName(), proc->getInstName());

    log().config().setProc(proc->getProcName(),proc->getInstName());
    proc->finish(param);
    log().config().disableProc();
  }
}

