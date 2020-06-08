// TreeManager.h
// Manager trees

#ifndef G2RECO_TREEMANAGER_H
#define G2RECO_TREEMANAGER_H

#include <typeinfo>
#include <typeindex> // require c++11...
#include <vector>

#include "TTree.h"
#include "TClass.h"
#include "TFile.h"

#include "Parameter.h"

namespace TreeProc{

  class TreeManager{

    class type_index2 : public std::type_index{
    public:
      type_index2(): type_index(typeid(int)){}
      type_index2(const std::type_info &id): type_index(id){}
    };

    private:
      TreeManager(); // singleton: create only by Instance()
    public:
      ~TreeManager(){}

    public:
      static TreeProc::TreeManager * instance();

      void config(Parameter *param){_param = param;}
      void initializeTree(); // automatic initialization with parameters; must call after config()
      void initializeTree(TTree *inputTree, const char *outputFileName, const char *outputTreeName = "Tree");
      
      // branches
      template <typename T> const T* getBranch(const char *name);
      template <typename T> const std::vector<const T *> & getBranchVec(const char *name);

      // register: common for obj and vec
      template <typename T> void registerBranch(const char *name, T * obj);
      
      // operations
      int getEntries()const {return _inTree->GetEntries();}
      bool getEntry(int i);
      void fill(){_outTree->Fill();}
      void write(){_outFile->cd();_outTree->Write();}
     
    private:
      static TreeProc::TreeManager *_instance;

      void * getBranch(const char *name, const std::type_info &type);
      
      // global maps
      typedef std::map<std::string, std::pair<type_index2, void *> > branchMap_type;
      branchMap_type _branchMap;
      typedef branchMap_type::iterator branchMap_iterator;

      Parameter *_param;
      // main tree
      TFile * _outFile;
      TTree * _outTree;
      TTree * _inTree;
      
      // previously read event number
      int _evIn;
  };

  // template implementation (inline)
  template <typename T> const T* TreeManager::getBranch(const char *name)
  {
    const T * pobj = static_cast<const T*> (getBranch(name, typeid(T)));
    if(pobj)
      return pobj;
    else
      throw("TreeManager::getBranch(): branch not found."); // message class should be developed
  }

  template <typename T> const std::vector<const T *> & TreeManager::getBranchVec(const char *name)
  {
    std::vector<const T *> * pvec = static_cast<std::vector<const T *> *>
      (getBranch(name, typeid(std::vector<const T*>)));
 
    if(pvec)
        return *pvec;
    else
      throw("TreeManager::getBranchVec(): branch not found."); // message class should be developed
  }

  template <typename T> void TreeManager::registerBranch(const char *name, T * pobj)
  {
    if(_branchMap.find(name) != _branchMap.end())
      throw("TreeManager::createBranch: the specified branch already exists.");
        
    _outTree->Branch(name, pobj);
    _branchMap[name] = std::pair<type_index2, void *>(typeid(T), pobj);
  }

/*
  template <typename T> T* TreeManager::createBranch(const char *name)
  {
    if(_branchMap.find(name))
      throw("TreeManager::createBranch: the specified branch already exists.");
        
    //T * pobj = static_cast<T *>(_branchMap[typeid(T)]->get());
    T * pobj = static_cast<T *>(TClass::GetClass(typeid(T))->New());
    _outTree->Branch(name, pobj);
    _branchMap[name] = std::pair<type_index2, void *>(typeid(T*), pobj);
  }
*/

}


#endif
