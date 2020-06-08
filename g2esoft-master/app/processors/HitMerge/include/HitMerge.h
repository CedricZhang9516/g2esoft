// HitMerge.h

#ifndef HitMerge_h
#define HitMerge_h

#include "TreeProc/Processor.h"
#include "app/objects.h"

#include <vector>
#include <string>

namespace g2esoft{
  class HitMerge : public TreeProc::Processor {
  public:
    HitMerge(const char *procName, const char *instName);
    virtual ~HitMerge();
    
    virtual void config(TreeProc::Parameter *);
    virtual void process(int nEvent, TreeProc::Parameter *);
    virtual void finish(TreeProc::Parameter *);

  private:
    std::vector<const StripCluster *> _stripClusters;
    
    std::string _parStripClusterName;
    std::string _parStripHitName;
    double _parTimeMatch;
  };
}

#endif
