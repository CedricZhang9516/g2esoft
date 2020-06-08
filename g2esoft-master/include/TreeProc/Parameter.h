// Parameter.h

#ifndef TREEPROC_PARAMETER_H
#define TREEPROC_PARAMETER_H

#include <map>
#include <string>
#include <typeindex>
#include <vector>

#include <iostream>

namespace TreeProc{

  class Parameter{
    friend class ParameterWriter;
    
    public:
      Parameter() : _parent(0){}
      ~Parameter(){}

    public:
      // access functions
      template <typename T>
        const T & get(const char *name, const T & defval = T())const{
          // checking existence
          paramType::const_iterator it = _param.find(name);
          if(it == _param.end()){
            // look for parent
            if(_parent){return _parent->get(name, defval);}
            else {
              return defval;
            }
          }
          // found: check type
          if(it->second.first != typeid(T)){
            throw("Parameter: type invalid.");
          }
          // convert to type T
          T * pret = static_cast<T*>(it->second.second);

          return *pret;
        }
      const char *name()const{return _name.c_str();}
      
      // non-const fucnctions: processors cannot access this
      Parameter * getChild(const char *name){
        //std::cout << "Parameter::getChild(): name = " << name << ", size = " << _children.size() << std::endl;
        
        auto it = _children.find(name);
        if(it != _children.end()) return it->second;
        else return 0;
      }
      Parameter * getParent(){return _parent;}

      typedef std::map<std::string, std::pair<std::type_index, void *> > paramType;
    private:
      paramType _param;
      
      std::string _name;
      
      Parameter * _parent;
      std::map<std::string, Parameter *> _children;
  };
  
  // class for right to write to parameter class
  // to be used for XMLParameterParser etc., not for general processors
  class ParameterWriter{
  protected:
    Parameter::paramType &getMap(Parameter *param){return param->_param;}
    void setParent(Parameter *param, Parameter *parent){param->_parent = parent;}
    void addChild(Parameter *param, Parameter *child){param->_children.insert(std::make_pair(child->name(), child));}
    void setName(Parameter *param, const char *name){param->_name = name;}
  };

  
}


#endif


