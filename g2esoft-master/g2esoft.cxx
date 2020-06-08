// g2esoft.cxx - main executable

// TreeProc includes
#include "TreeProc/ProcManager.h"
#include "TreeProc/TreeManager.h"
#include "TreeProc/RunManager.h"
#include "TreeProc/LogStream.h"

// test purpose
#include "TreeProc/XMLParser.h"
#include <vector>
#include <string>

// ROOT includes
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TClassTable.h"

// STL includes
#include <iostream>
#include <exception>

using namespace std;
using namespace TreeProc;

int main(int argc, char *argv[])
{
  try{
    // get xml file before initializing TApplication since it takes out xml filename!
    string xmlfile = "test.xml";
    if(argc > 1)
      xmlfile = argv[1];

    TApplication app("g2esoft",&argc, argv);
  
    XMLParser parser;
    vector<string> libs;
    vector<pair<string, string> > procs;
    Parameter param;
    
    cout << "Reading steering file: " << xmlfile << endl;
    log() << "Reading steering file: " << xmlfile << endl;
    
    parser.loadFile(xmlfile.c_str());
    parser.getLibraries(libs);
    parser.getProcessors(procs);
    parser.getParameters(&param);

    // load library: need LD_LIBRARY_PATH
    RunManager::instance()->loadLibraries(libs);

    //gClassTable->Print();

    // initialization of logging
    log().config().initializeFromParameter(&param);
    
    // initialization of trees
    TreeManager::instance()->config(&param);
    TreeManager::instance()->initializeTree();

    // initialization of processors
    ProcManager::instance()->initializeProcs(procs);
    ProcManager::instance()->config(&param);

    RunManager::instance()->config(&param, parser.getArranger());
    RunManager::instance()->eventLoop();

  }catch(exception &e){
    // if log is still temporal string buffer output the contents to cerr
    log().config().openFile("cerr",true);
    cout << "Exception occurred: " << e.what() << endl;
  }
  
  return 0;
}
