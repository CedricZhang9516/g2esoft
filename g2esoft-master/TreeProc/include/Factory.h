// Factory.h
// Factory & FactoryBase class
// Please note that the 

#ifndef TREEPROC_FACTORY_H
#define TREEPROC_FACTORY_H

#include <map>
#include <string>
#include <string.h>

namespace TreeProc{

  template<class T> class ObjectFactory;
  
  template <class T> class FactoryBase{
    protected:
      FactoryBase(){}
    public:
      virtual T * get(const char *name) = 0; // pure virtual OK?
      virtual ~FactoryBase(){}
  };
  
  template<class base, class T> class Factory : public FactoryBase<base>{
    public:
      Factory(const char *name){TreeProc::ObjectFactory<base>::instance()->registerFactory(name,this);_procName = name;}
      ~Factory(){}
      virtual base * get(const char *name) {return static_cast<base *>(new T(_procName.c_str(), name));}
      
    private:
      std::string _procName;
  };
  
  template<class T> class ObjectFactory{
    typedef typename std::map<std::string, TreeProc::FactoryBase<T> *> FactoryMapType;

  private:
    ObjectFactory(){}
  public:
    ~ObjectFactory(){}

    static ObjectFactory<T> * instance(){if(_instance)return _instance; else {_instance = new ObjectFactory<T>; return _instance;}}
    void registerFactory(const char *name, TreeProc::FactoryBase<T> *factory){
      if(_factoryMap.find(name) != _factoryMap.end()){
        throw("ObjectFactgory::registerFactory(): factory already registered.");
      }
      _factoryMap[name] = factory;
    }
    
    T * makeObject(const char *procName, const char *instanceName){
      //FactoryMapType::iterator it;
      typename FactoryMapType::iterator it;
      if((it = _factoryMap.find(procName)) == _factoryMap.end()){
        throw("ObjectFactory::makeObject(): object cannot found in the factory.");
      }
      if(strlen(instanceName) == 0)
        instanceName = procName;
  
      return it->second->get(instanceName);
    }
    
    
  private:
    static ObjectFactory<T> * _instance;
   
    FactoryMapType _factoryMap;
    
  };
  
  template<class T> ObjectFactory<T> * ObjectFactory<T>::_instance = 0;
}


#endif
