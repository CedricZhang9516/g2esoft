// TestProc.cxx

#include "TreeProc/TreeManager.h"
#include "TreeProc/LogStream.h"
#include "TestProc.h"

#include <iostream>

using namespace g2esoft;
using namespace TreeProc;
using namespace std;

static Factory<Processor, TestProc> _f("TestProc");

TestProc::TestProc(const char *procName, const char *instName) : Processor(procName, instName)
{}

TestProc::~TestProc(){}

void TestProc::config(Parameter *param){
  _parInputFileName = param->get("InputFile", string("InputFile_default"));
  _parMCParticleName = param->get("MCParticleName",string("MCParticles"));
  _parTestString = param->get("test", string("Test_default"));
  _parTestStrings = param->get("test2", vector<string>{string("Test1"),string("Test2")});
}

void TestProc::process(int nEvent, Parameter *){
  log() << "Come to event " << nEvent << endl;
  
  const vector<const MCParticle *> &mc = TreeManager::instance()->getBranchVec<MCParticle>(_parMCParticleName.c_str());
  log() << "Number of MCP = " << mc.size() << endl;
  
  // parameter test
  log() << "Parameter test: " << _parTestString << endl;
  
  log() << "Parameter test2: " << endl;
  for(unsigned int i=0;i<_parTestStrings.size();i++)
    log() << _parTestStrings.at(i) << endl;
  log() << "Parameter test2 end." << endl;
  log() << "Parameter InputFile: " << _parInputFileName << endl;
  
}
void TestProc::finish(Parameter *){}

