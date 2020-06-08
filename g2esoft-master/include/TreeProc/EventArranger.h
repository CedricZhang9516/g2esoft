// EventArranger.h
// The base class to inherite to make own processors

#ifndef TREEPROC_EVENTARRANGER_H
#define TREEPROC_EVENTARRANGER_H

#include "TreeProc/ProcManager.h" // necessary for template functions
#include "TreeProc/Parameter.h"

#include <string>

namespace TreeProc{

  class EventArranger{
    protected:
      EventArranger(const char *procName, const char *instName); // create only by Factory::get()
    public:
      virtual ~EventArranger(){}

    public:
      // access functions
      virtual void config(Parameter *){}
      virtual bool arrange(int nEventIn, Parameter *) = 0; // true: arrange succeeded (to process) false: arrange failed (to continue)
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
