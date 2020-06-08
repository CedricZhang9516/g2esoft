// ProcManager.h
// Manage processors

#ifndef G2RECO_PROCMANAGER_H
#define G2RECO_PROCMANAGER_H

#include "Factory.h"
#include <vector>
#include <string>
#include "Parameter.h"

namespace TreeProc{
  // forward declarations
  class Processor;
  class Parameter;
  
  class ProcManager{
    private:
      ProcManager(); // singleton: create only by Instance()
    public:
      ~ProcManager(){}

    public:
      static TreeProc::ProcManager * instance();

      void initializeProc(const char *procName, const char *instanceName = "");
      void initializeProcs(std::vector<std::pair<std::string, std::string> > &procs){
        for(unsigned int i=0;i<procs.size();i++){
          initializeProc(procs[i].first.c_str(), procs[i].second.c_str());
        }
      }

      void config(Parameter *rootParams);
      void process(int nEvent);
      void finish();
      
      inline Parameter * getParam(const std::string &procName, const std::string &instName);

  private:
      static TreeProc::ProcManager *_instance;

      
      // global maps
      typedef std::pair<std::string, TreeProc::Processor *> procPair;
      typedef std::vector< procPair > procVecType;
      typedef std::map<std::string, TreeProc::Processor *> procMapType;
      // duplicate info: but keep map to find from name
      procVecType _processors;
      procMapType _procMap;
      
      Parameter *_rootParams;
  };
  
  inline Parameter * ProcManager::getParam(const std::string &procName, const std::string &instName)
  {
    Parameter *parProc = _rootParams->getChild(procName.c_str());
    if(parProc == 0)return _rootParams;
    
    Parameter *parInst = parProc->getChild(instName.c_str());
    if(parInst == 0)return parProc;
    return parInst;
  }

}

#endif
