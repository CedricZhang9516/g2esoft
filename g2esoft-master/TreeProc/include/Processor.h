// Processor.h
// The base class to inherite to make own processors

#ifndef G2RECO_PROCESSOR_H
#define G2RECO_PROCESSOR_H

#include "TreeProc/ProcManager.h" // necessary for template functions
#include "TreeProc/Parameter.h"

#include <string>

namespace TreeProc{

  class Processor{
    protected:
      Processor(const char *procName, const char *instName); // create only by Factory::get()
    public:
      virtual ~Processor(){}

    public:
      // access functions
      virtual void config(Parameter *){}
      virtual void process(int /*nEvent*/, Parameter *){}
      virtual void finish(Parameter *){}
      
      // accessor to names
      const char * getProcName()const{return _procName.c_str();}
      const char * getInstName()const{return _instName.c_str();}

    private:
      std::string _instName; // instance name
      std::string _procName; // processor name
  };  
}


#endif
