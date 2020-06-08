// RunManager.h
// Event loop holder

#ifndef G2RECO_RUNMANAGER_H
#define G2RECO_RUNMANAGER_H

#include <vector>
#include <string>

#include "Parameter.h"
#include "EventArranger.h"

namespace TreeProc{

  class RunManager{
    private:
      RunManager(); // singleton: create only by Instance()
    public:
      ~RunManager(){}

    public:
      static TreeProc::RunManager * instance();

      void loadLibraries(std::vector<std::string> &libs);
      void config(Parameter *param, const char *arranger);
      void eventLoop();
      int getNextEntry();

    private:
      static TreeProc::RunManager *_instance;

      const Parameter *_param;
      EventArranger * _arranger;
      int _evIn;
      int _evOut;
      
      // steering parmeters
      int _parEventNumToProcess; // event number to be processed
      int _parEventIdxToStart; // event index to start
      int _parWritePeriod; // number of output event to accumurate before writing to tree
      int _parPrintPeriod; // number of input event before notice in log output
      const char * _parArrangerName; // name of arranger: not at global parameter but at run tag
      std::string _parInputFile; // input file name
  };
}

#endif
