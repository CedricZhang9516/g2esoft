// TreeManager.cxx

// ROOT includes
#include "TFile.h"

// TreeProc includes
#include "TreeProc/TreeManager.h"
#include "TreeProc/LogStream.h"

#include <iostream>
#include <exception>
#include <stdexcept>

using namespace std;
using namespace TreeProc;

// singleton
TreeManager * TreeManager::_instance = 0;
TreeManager * TreeManager::instance(){
  if(_instance) return _instance;
  else return _instance = new TreeManager;
}

TreeManager::TreeManager()
: _inTree(0)
, _outTree(0)
, _outFile(0)
, _evIn(-1)
{}

void TreeManager::initializeTree(TTree *inputTree, const char *outputFileName, const char *outputTreeName)
{
  _inTree = inputTree;
  _outFile = TFile::Open(outputFileName,"new");
  if(!_outFile){
    log("fatal") << "TreeManager::initializeTree: Output file " << outputFileName << " open failed. check if the file already exists." << endl;
    throw runtime_error("Output file open failed.");
  }
  
  const char *treeName = outputTreeName;
  if(outputTreeName == "") treeName = _inTree->GetName();
  
  _outTree = new TTree(treeName, "");
  _outTree->BranchRef();
  if(_inTree)
    _outTree->AddFriend(_inTree);
}

void TreeManager::initializeTree()
{
  string inputFile = _param->get("InputFile", string("testIn.root"));
  log() << "Input file: " << inputFile << endl;

  /*
  TFile *file = TFile::Open(inputFile.c_str());
  TTree *tree = 0;
  if(!file){
    log() << "File " << inputFile << " not found!" << endl;
  }else{
    string inputTree = _param->get("InputTree", string("trk"));
    tree = static_cast<TTree*>(file->Get(inputTree.c_str()));
    if(!tree){
      log() << "Tree " << inputTree << " is not found in the file " << inputFile << endl;
    }
  }
  */
  TTree *tree = 0;

  if(inputFile != ""){
    TFile *file = TFile::Open(inputFile.c_str());
    if(!file){
      log() << "File " << inputFile << " not found!" << endl;
    }else{
      string inputTree = _param->get("InputTree", string("trk"));
      tree = static_cast<TTree*>(file->Get(inputTree.c_str()));
      if(!tree){
        log() << "Tree " << inputTree << " is not found in the file " << inputFile << endl;
      }
    }
  }else{
    tree = new TTree ("trk", "trk");
  }
  
  string outputFile = _param->get("OutputFile", string("testOut.root"));
  log() << "Output file: " << outputFile << endl;
  string outputTree = _param->get("OutputTree", string("trk"));
  TreeManager::instance()->initializeTree(tree,outputFile.c_str(),outputTree.c_str());
}

void * TreeManager::getBranch(const char *name, const type_info &type)
{
  branchMap_iterator it;
  if((it = _branchMap.find(name)) != _branchMap.end()){
    // type confirmation
    if(it->second.first != type){
      string message = "TreeManager::getBranch(): type of the branch name not correct for ";
      message += name;
      message += ".";
      throw runtime_error(message);
      //throw runtime_error("TreeManager::getBranch(): type of the branch name not correct.");
    }
    return it->second.second;
  }
  else{
    // not found in the existing map: try SetBranchAddress
    TBranch *br = _outTree->FindBranch(name);
    if(!br){
      string message = "TreeManager::getBranch(): branch not found for ";
      message += name;
      message += ".";
      throw runtime_error(message);
      //throw runtime_error("TreeManager::getBranch(): branch not found.");
    }
    // create object
    log() << "Create object for " << TClass::GetClass(type)->GetName() << endl;
    void *ptr = TClass::GetClass(type)->New();
    _branchMap[name] = std::pair<type_index2, void *>(type, ptr);
    // be careful! should be void *& (address is regietered in tree)
    void *&ptr2 = _branchMap[name].second;
    
    // SetBranchAddress
    _outTree->SetBranchAddress(name, &ptr2, TClass::GetClass(type), kOther_t, true);
    if(_evIn >= 0)br->GetEntry(_evIn);
    return _branchMap[name].second;
  }
}

bool TreeManager::getEntry(int i)
{
  if(_inTree->GetEntry(i) < 0){
    log("error") << "getEntry(): Failed to move to new event " << i << "." << endl;
    return false;
  }
  _evIn = i;
  
  return true;
}
