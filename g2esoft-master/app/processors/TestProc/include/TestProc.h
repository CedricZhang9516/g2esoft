// TestProc.h

#ifndef TestProc_h
#define TestProc_h

#include "TreeProc/Processor.h"
#include "app/objects.h"

namespace g2esoft{
  class TestProc : public TreeProc::Processor {
  public:
    TestProc(const char *procName, const char *instName);
    virtual ~TestProc();
    
    virtual void config(TreeProc::Parameter *);
    virtual void process(int nEvent, TreeProc::Parameter *);
    virtual void finish(TreeProc::Parameter *);

  private:
    std::string _parInputFileName;
    std::string _parMCParticleName;
    std::string _parTestString;
    std::vector<std::string> _parTestStrings;
  };
}

#endif
