// XMLParser.cxx

#include "TreeProc/XMLParser.h"
#include "TreeProc/LogStream.h"

#include <string.h>

#include <sstream>
#include <iostream>

#include <TXMLAttr.h>
#include <TList.h>

using namespace std;
using namespace TreeProc;

void XMLParser::loadFile(const char *file)
{
  _parser = new TDOMParser;
  _parser->SetValidate(false);
  int ret = _parser->ParseFile(file);
  
  if(ret < 0)throw("XML file readout failed.");
  
  _rootNode = _parser->GetXMLDocument()->GetRootNode();
  
}

void XMLParser::getLibraries(std::vector<std::string> &libs)
{
  // look for "libraries" tag
  TXMLNode *node = _rootNode->GetChildren();
  
  log() << "Libraries:" << endl;
  while(node){
    if(!strcasecmp("libraries", node->GetNodeName())){
      istringstream is(string(node->GetText()));
      string s;
      while(1){
        s.clear();
        is >> s;
        if(s.empty())break;
        log() << s << endl;
        libs.push_back(s);
      }
    }
    node = node->GetNextNode();
  }
  log() << "Libraries end." << endl << endl;
}

const char *XMLParser::getArranger()
{
  TXMLNode *node = _rootNode->GetChildren();
  while(node){
    if(!strcasecmp("run", node->GetNodeName())){
      TList *attrList = node->GetAttributes();
      if(!attrList){return 0;} // no arranger supplied in run parameter
      TXMLAttr *attr = static_cast<TXMLAttr *>(attrList->FindObject("arranger"));
      if(!attr){return 0;} // no arranger element
      return attr->GetValue();
    }
    node = node->GetNextNode();
  }
  return 0;
}

void XMLParser::getProcessors(std::vector<std::pair<std::string, std::string> > &procs)
{
  log() << "Run:" << endl;

  // look for "run" tag
  TXMLNode *node = _rootNode->GetChildren();
  while(node){
    if(!strcasecmp("run", node->GetNodeName())){
      TXMLNode *cnode = node->GetChildren();
      while(cnode){
        if(!strcasecmp("processor", cnode->GetNodeName())){
          string procname, instname;
          
          TList *attrList = cnode->GetAttributes();
          if(!attrList){throw("XMLParser: processor element has no name.");}
          TXMLAttr *attr = static_cast<TXMLAttr *>(attrList->FindObject("name"));
          if(!attr){throw("XMLParser: processor element has no name.");}
          procname = attr->GetValue();
          
          attr = static_cast<TXMLAttr *>(attrList->FindObject("instance"));
          if(attr)instname = attr->GetValue();
          else instname = procname;

          log() << "Proc name: " << procname << ", inst name: " << instname << endl;
          procs.push_back(make_pair(procname, instname));
        }
        cnode = cnode->GetNextNode();
      }
    }
    node = node->GetNextNode();
  }
  log() << "Run end." << endl << endl;
}

void XMLParser::getParameters(Parameter *rootParam)
{
  log() << "Parameters:" << endl;

  // look for "parameters" tag
  TXMLNode *node = _rootNode->GetChildren();
  while(node){
    if(!strcasecmp("parameters", node->GetNodeName())){
      TXMLNode *cnode = node->GetChildren();
      getParamElement(rootParam, cnode);
    }
    node = node->GetNextNode();
  }
  log() << "Parameters end." << endl << endl;
}

void XMLParser::getParamElement(Parameter *param, TXMLNode *node)
{
  while(node){
    if(!strcasecmp("bool", node->GetNodeName())){
      getParamNode<bool>(param, node);
    }
    else if(!strcasecmp("int", node->GetNodeName())){
      getParamNode<int>(param, node);
    }
    else if(!strcasecmp("double", node->GetNodeName())){
      getParamNode<double>(param, node);
    }
    else if(!strcasecmp("string", node->GetNodeName())){
      getParamNode<string>(param, node);
    }
    else if(!strcasecmp("boolVec", node->GetNodeName())){
      getParamNode<vector<bool> >(param, node);
    }
    else if(!strcasecmp("intVec", node->GetNodeName())){
      getParamNode<vector<int> >(param, node);
    }
    else if(!strcasecmp("doubleVec", node->GetNodeName())){
      getParamNode<vector<double> >(param, node);
    }
    else if(!strcasecmp("stringVec", node->GetNodeName())){
      getParamNode<vector<string> >(param, node);
    }
    else if(!strcasecmp("processor", node->GetNodeName())
      || !strcasecmp("arranger", node->GetNodeName())
      || !strcasecmp("instance", node->GetNodeName()) ){
      // make a child parameter
      Parameter *child = new Parameter;
      // set parameter name
      TList *attrList = node->GetAttributes();
      if(!attrList){throw("XMLParser: parameter element has no name.");}
      TXMLAttr *attr = static_cast<TXMLAttr *>(attrList->FindObject("name"));
      if(!attr){throw("XMLParser: parameter element has no name.");}
      setName(child, attr->GetValue());

      setParent(child, param);
      addChild(param, child);
    
      log() << "New parameter added: name =  " << attr->GetValue() << endl;
      
      // recuvsive call for child parameter
      getParamElement(child, node->GetChildren());

      log() << "Parameter " << attr->GetValue() << " read. Number of parameters = " << getMap(child).size() << endl;
    }
    
    node = node->GetNextNode();
  }
}

// conversion of non-vectors
template <typename T> T * XMLParser::Converter<T>::operator() (const char *text){
  istringstream is; // I don't know why, but if I write is(string(text)) the compile fails at operator>> while it's OK with non-template usage.
  is.str(text);
  
  T * val = new T;
  is >> *val;
  
  return val;
}

// specialization for string; just assign
string * XMLParser::Converter<string>::operator() (const char *text){
  string *val = new string(text);
  return val;
}

// conversion of vectors
template <typename T> vector<T> * XMLParser::Converter<vector<T> >::operator() (const char *text){
  istringstream is;
  is.str(text);
  
  vector<T> *vals = new vector<T>;

  while(is.good()){
    T val;
    is >> val;
    vals->push_back(val);
  }
  
  return vals;
}

template <typename T> void XMLParser::getParamNode(Parameter *param, TXMLNode *node)
{
  // look for name
  TList *attrList = node->GetAttributes();
  if(!attrList){throw("XMLParser: parameter element has no name.");}
  TXMLAttr *attr = static_cast<TXMLAttr *>(attrList->FindObject("name"));
  if(!attr){throw("XMLParser: parameter element has no name.");}
  
  string paramName = attr->GetValue();

  T *val;
  const char *text = node->GetText();
  if(text != 0)
    // need a conversion from const char * to T
    val = Converter<T>()(text);
  else
    val = Converter<T>()("");

  // insert to parameter class
  Parameter::paramType &map = getMap(param);
  map.insert(make_pair(paramName, make_pair<type_index, void *>(typeid(T), val)));
}

